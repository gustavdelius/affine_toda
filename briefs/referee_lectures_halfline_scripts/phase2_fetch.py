from pathlib import Path
import requests,tarfile,io,gzip,concurrent.futures
base=Path('referee_scratch/lectures_halfline/sources')
ids='9404108 9407148 9501098 9408004 9306089 9311027 9304141 9506031 9609233 9707235 9709078 9806003 9807189 9809140 0106070 0112023 0507053 0906.1215'.split()
def fetch(i):
 d=base/i
 if d.exists(): return i+' exists'
 aid=i if '.' in i else 'hep-th/'+i
 try:
  r=requests.get('https://arxiv.org/src/'+aid,timeout=45);r.raise_for_status()
  d.mkdir()
  try:
   with tarfile.open(fileobj=io.BytesIO(r.content),mode='r:*') as t:t.extractall(d,filter='data')
  except tarfile.ReadError:
   try:b=gzip.decompress(r.content)
   except OSError:b=r.content
   (d/'source.tex').write_bytes(b)
  return i+' OK'
 except Exception as e:return i+' ERROR '+str(e)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
 for s in ex.map(fetch,ids):print(s,flush=True)
