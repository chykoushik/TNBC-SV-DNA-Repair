import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
LIMIT = 1_000_000_000

def safe_path(root, relative):
    root = Path(root).resolve()
    relative = Path(relative)
    if relative.is_absolute() or '..' in relative.parts:
        raise ValueError('Unsafe relative path')
    result = (root / relative).resolve()
    if not result.is_relative_to(root) or result == root:
        raise ValueError('Path escapes permitted root')
    return result

def hashes(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            sha.update(block); md5.update(block)
    return sha.hexdigest(), md5.hexdigest()

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def record(root, entry):
    p = safe_path(root, 'outputs/acquisition/download_manifest.jsonl')
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open('a', encoding='utf-8') as f:
        f.write(json.dumps(entry, sort_keys=True) + '\n')

def expand(root, resource):
    if resource.get('snapshot'):
        path = safe_path(root, resource['snapshot'])
        if hashes(path)[0] != resource['snapshot_sha256']:
            raise ValueError('GDC inventory changed; explicitly repin before acquisition')
        entries = json.loads(path.read_text(encoding='utf-8'))['files']
        result = []
        for x in entries:
            if x['access'] != 'open':
                raise ValueError('Controlled file rejected')
            name = x['file_name']
            if '/' in name or '\\' in name or ':' in name:
                raise ValueError('Unsafe remote filename')
            result.append(dict(resource, url='https://api.gdc.cancer.gov/data/' + x['file_id'],
                               filename=x['file_id'] + '/' + name,
                               expected_bytes=x['file_size'], md5=x['md5sum']))
        return result
    return [resource]

def transfer(root, item, *, resume=False, overwrite=False, opener=urllib.request.urlopen):
    if not item['url'].startswith('https://'):
        raise ValueError('Only HTTPS public source URLs allowed')
    expected = item.get('expected_bytes')
    if not isinstance(expected, int) or expected <= 0:
        raise ValueError('Unknown or invalid size: manual acquisition required')
    external = safe_path(root, 'data/external')
    destination = safe_path(external, item['id'] + '/' + item['filename'])
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + '.part')
    state = destination.with_name(destination.name + '.part.json')
    lock = destination.with_name(destination.name + '.lock')
    for p in (partial, state, lock):
        if p.is_symlink() or not p.resolve().is_relative_to(external):
            raise ValueError('Unsafe sidecar path')
    entry = dict(dataset_id=item['id'], path=destination.relative_to(Path(root).resolve()).as_posix(),
                 source_url=item['url'], provider=item['provider'], release=item['release'],
                 started_utc=now(), expected_bytes=expected, provider_md5=item.get('md5'))
    lock_handle = None
    try:
        lock_handle = lock.open('x')
        if destination.exists() and not overwrite:
            print(f"SKIP existing: {destination}")
            return {"status": "skipped", "path": str(destination)}
        offset = 0
        old = {}
        if partial.exists() and not overwrite:
            if not resume or not state.exists():
                raise FileExistsError('Partial exists; use --resume or explicit --overwrite')
            old = json.loads(state.read_text(encoding='utf-8'))
            if old.get('url') != item['url'] or old.get('expected_bytes') != expected:
                raise ValueError('Partial provenance mismatch')
            offset = partial.stat().st_size
            if not 0 < offset < expected:
                raise ValueError('Partial size cannot be resumed; use explicit --overwrite')
            if not old.get('validator'):
                raise ValueError('No ETag/Last-Modified validator; safe resume unavailable')
        headers = {'User-Agent': 'TNBC-V2-public-acquisition/1.0', 'Accept-Encoding': 'identity'}
        if offset:
            headers.update({'Range': f'bytes={offset}-', 'If-Range': old['validator']})
        request = urllib.request.Request(item['url'], headers=headers)
        with opener(request, timeout=120) as response:
            status = response.status
            if not response.geturl().startswith('https://'):
                raise ValueError('Insecure redirect rejected')
            if offset:
                match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', response.headers.get('Content-Range', ''))
                if status != 206 or not match or tuple(map(int, match.groups())) != (offset, expected-1, expected):
                    raise ValueError('Server did not honor exact resume range; partial preserved')
            elif status != 200:
                raise ValueError(f'Unexpected HTTP status {status}')
            length = response.headers.get('Content-Length')
            if length is not None and int(length) != expected-offset:
                raise ValueError('Remote size differs from pinned size')
            content_type = response.headers.get('Content-Type', '').lower()
            if 'text/html' in content_type:
                raise ValueError('HTML verification/error page rejected')
            etag = response.headers.get('ETag')
            validator = etag if etag and not etag.startswith('W/') else response.headers.get('Last-Modified')
            if offset and validator and old['validator'] != validator:
                raise ValueError('Remote validator changed')
            state.write_text(json.dumps(dict(url=item['url'], expected_bytes=expected, validator=validator)), encoding='utf-8')
            count = offset
            with partial.open('ab' if offset else 'wb') as f:
                while True:
                    block = response.read(1024 * 1024)
                    if not block: break
                    if count == 0 and block.lstrip().lower().startswith((b'<!doctype html', b'<html')):
                        raise ValueError('HTML payload rejected')
                    count += len(block)
                    if count > expected:
                        raise ValueError('Transfer exceeded pinned size')
                    f.write(block)
            entry['final_url'] = response.geturl()
        if count != expected:
            raise ValueError(f'Truncated transfer: {count} != {expected}')
        sha, md5 = hashes(partial)
        if item.get('md5') and md5.lower() != item['md5'].lower():
            raise ValueError('Provider MD5 mismatch')
        if destination.exists() and not overwrite:
            raise FileExistsError('Destination appeared during transfer')
        if overwrite:
            os.replace(partial, destination)
        else:
            # Exclusive destination creation
            # Atomic verified publication
            os.link(partial, destination)
            partial.unlink()
        state.unlink(missing_ok=True)
        entry.update(status='success', size_bytes=count, sha256=sha, md5=md5, completed_utc=now())
        record(root, entry)
        return entry
    except Exception as exc:
        entry.update(status='failed', error=str(exc), completed_utc=now())
        record(root, entry)
        raise
    finally:
        if lock_handle is not None:
            lock_handle.close(); lock.unlink(missing_ok=True)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', action='append', default=[])
    parser.add_argument('--recommended', action='store_true')
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--allow-large', action='store_true', help='Explicitly authorize selected total >=1 GB')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    config = json.loads((ROOT / 'config/data_sources.yaml').read_text(encoding='utf-8'))
    available = {r['id']: r for r in config['resources']}
    unknown = set(args.dataset) - set(available)
    if unknown: parser.error('Unknown, manual, local-only or controlled IDs: ' + ', '.join(sorted(unknown)))
    selected = [r for r in available.values() if r['id'] in args.dataset or (args.recommended and r['recommended'])]
    if not selected:
        if args.download: parser.error('Select --dataset or --recommended explicitly')
        selected = list(available.values())
    files = [f for r in selected for f in expand(ROOT, r)]
    total = sum(f['expected_bytes'] for f in files)
    print(json.dumps(dict(mode='download' if args.download else 'dry_run', datasets=[r['id'] for r in selected], files=len(files), bytes=total), indent=2))
    if not args.download: return 0
    if total >= LIMIT and not args.allow_large:
        parser.error('Selected transfer totals >=1 GB; requires --allow-large. Nothing downloaded.')
    for i, f in enumerate(files, 1):
        print(f"[{i}/{len(files)}] {f['id']}: {f['filename']}", flush=True)
        transfer(ROOT, f, resume=args.resume, overwrite=args.overwrite)
    return 0

if __name__ == '__main__':
    try: sys.exit(main())
    except Exception as exc:
        print(f'Acquisition stopped: {exc}', file=sys.stderr); sys.exit(1)
