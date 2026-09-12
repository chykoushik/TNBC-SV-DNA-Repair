import json
import sys
from download_public_data import ROOT, hashes, safe_path, now

def verify(root):
    manifest = safe_path(root, 'outputs/acquisition/download_manifest.jsonl')
    if not manifest.exists():
        return {'checked_utc': now(), 'status': 'no_download_manifest', 'files': [], 'passed': False}
    latest = {}
    for line in manifest.read_text(encoding='utf-8').splitlines():
        item = json.loads(line)
        if item.get('status') == 'success': latest[item['path']] = item
    results = []
    for relative, expected in latest.items():
        try:
            path = safe_path(root, relative)
            if not path.is_relative_to(safe_path(root, 'data/external')):
                raise ValueError('Manifest payload is outside data/external')
            sha, md5 = hashes(path)
            good = path.stat().st_size == expected['size_bytes'] and sha == expected['sha256']
            if expected.get('provider_md5'): good = good and md5 == expected['provider_md5'].lower()
            results.append(dict(path=relative, passed=good, sha256=sha))
        except Exception as exc:
            results.append(dict(path=relative, passed=False, error=str(exc)))
    return dict(checked_utc=now(), status='verified' if latest else 'no_successful_downloads', files=results, passed=bool(results) and all(r['passed'] for r in results))

if __name__ == '__main__':
    report = verify(ROOT)
    output = safe_path(ROOT, 'outputs/acquisition/download_verification.json')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))
    sys.exit(0 if report['passed'] else 1)
