import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from download_public_data import transfer, safe_path
import download_public_data as downloader
from verify_downloads import verify

class Response(io.BytesIO):
    def __init__(self, data, status=200, **headers):
        super().__init__(data); self.status=status
        self.headers={'Content-Length':str(len(data)), 'ETag':'"version1"', **headers}
    def geturl(self): return 'https://example.org/data'

class AcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.data=b'gene,sample\nA,1\n'
        self.item=dict(id='synthetic',filename='data.csv',url='https://example.org/data',expected_bytes=len(self.data),md5=hashlib.md5(self.data).hexdigest(),provider='synthetic test',release='1')
    def call(self, response=None, **kwargs):
        return transfer(self.root,self.item,opener=lambda *a,**k:response or Response(self.data),**kwargs)
    def partial(self):
        p=self.root/'data/external/synthetic/data.csv.part';p.parent.mkdir(parents=True)
        p.write_bytes(self.data[:5])
        p.with_name('data.csv.part.json').write_text(json.dumps(dict(url=self.item['url'],expected_bytes=len(self.data),validator='"version1"')))
        return p
    def test_success_and_tampering(self):
        self.call();self.assertTrue(verify(self.root)['passed'])
        (self.root/'data/external/synthetic/data.csv').write_bytes(b'corrupt')
        self.assertFalse(verify(self.root)['passed'])
    def test_never_overwrite(self):
        self.call()
        with self.assertRaises(FileExistsError):self.call()
        self.assertEqual((self.root/'data/external/synthetic/data.csv').read_bytes(),self.data)
    def test_resume(self):
        self.partial()
        self.call(Response(self.data[5:],206,**{'Content-Range':f'bytes 5-{len(self.data)-1}/{len(self.data)}'}),resume=True)
        self.assertTrue(verify(self.root)['passed'])
    def test_resume_ignored_preserves_partial(self):
        p=self.partial()
        with self.assertRaises(ValueError):self.call(resume=True)
        self.assertEqual(p.read_bytes(),self.data[:5])
    def test_wrong_range(self):
        self.partial()
        with self.assertRaises(ValueError):self.call(Response(self.data[5:],206,**{'Content-Range':f'bytes 4-{len(self.data)-1}/{len(self.data)}'}),resume=True)
    def test_truncated(self):
        with self.assertRaises(ValueError):self.call(Response(self.data[:-1],**{'Content-Length':str(len(self.data))}))
        self.assertFalse((self.root/'data/external/synthetic/data.csv').exists())
    def test_bad_checksum(self):
        self.item['md5']='0'*32
        with self.assertRaises(ValueError):self.call()
    def test_html(self):
        with self.assertRaises(ValueError):self.call(Response(self.data,**{'Content-Type':'text/html'}))
    def test_path_escape(self):
        with self.assertRaises(ValueError):safe_path(self.root,'../outside')
    def test_no_manifest(self):self.assertFalse(verify(self.root)['passed'])
    def test_explicit_overwrite(self):
        self.call();self.call(overwrite=True)
        self.assertTrue(verify(self.root)['passed'])
    def test_large_gate_before_network(self):
        with patch.object(sys,'argv',['download_public_data.py','--recommended','--download']), patch.object(downloader,'transfer') as network:
            with self.assertRaises(SystemExit) as result:downloader.main()
            self.assertEqual(result.exception.code,2);network.assert_not_called()
    def test_unknown_size(self):
        self.item['expected_bytes']=None
        with self.assertRaises(ValueError):self.call()
    def test_changed_validator(self):
        self.partial()
        with self.assertRaises(ValueError):self.call(Response(self.data[5:],206,**{'Content-Range':f'bytes 5-{len(self.data)-1}/{len(self.data)}','ETag':'"changed"'}),resume=True)

if __name__=='__main__':unittest.main()
