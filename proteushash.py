import hashlib
import base64
import sys

#http://hashdist.readthedocs.org/en/latest/core/source_cache.html#source-keys
hasher = hashlib.sha256()
f = open(sys.argv[1])
while True:
  chunk = f.read(1024)
  if not chunk:
    break
  hasher.update(chunk)
#http://hashdist.readthedocs.org/en/latest/core/hasher.html#hashdist.core.hasher.format_digest
key = base64.b32encode(hasher.digest()[:20]).lower()
print key
