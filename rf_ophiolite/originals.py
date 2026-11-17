"""Prove that the original data files never change.

'make' records a SHA-256 checksum of every file in the original data folders before
any processing; 'verify' checks every file against it later. A SHA-256 match means
the file is byte-for-byte identical: every sample and every header value.

    python rfo.py originals make   data/raw_sac data/raw_mseed -m originals_sha256.txt
    python rfo.py originals verify -m originals_sha256.txt

'make' also writes <manifest>.folders (the folder list, used to find new files).
'verify' reports changed, missing and new files and ends with
"RESULT: the original data are unchanged" (exit code 0) or a list of differences
(exit code 1).

The manifest is a standard checksum list, so it can also be checked without Python,
from the folder that holds it:
    shasum -a 256 -c --quiet originals_sha256.txt      (macOS)
    sha256sum -c --quiet originals_sha256.txt          (Linux)
This prints nothing if every listed file matches; it does not look for new files.

Hidden system files (names starting with '.') are ignored. Standard library only.
"""

import argparse
import hashlib
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

CHUNK = 8 * 1024 * 1024


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:  # read-only
        for block in iter(lambda: fh.read(CHUNK), b''):
            h.update(block)
    return h.hexdigest()


def same_or_inside(path, folder):
    """True if path is folder or below it (device/inode comparison: links and letter case cannot hide it)."""
    fst = os.stat(folder)
    q = os.path.abspath(path)
    while True:
        if os.path.exists(q) and os.path.samestat(os.stat(q), fst):
            return True
        parent = os.path.dirname(q)
        if parent == q:
            return False
        q = parent


def list_files(root):
    out = []
    for d, sub, files in os.walk(root, followlinks=True):
        sub[:] = sorted(s for s in sub if not s.startswith('.'))
        out += [os.path.join(d, f) for f in sorted(files) if not f.startswith('.')]
    return out


def hash_all(paths, workers):
    res, t0 = {}, time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for k, (p, h) in enumerate(zip(paths, ex.map(sha256, paths)), 1):
            res[p] = h
            if k % 5000 == 0 or k == len(paths):
                print('   %d / %d files hashed (%.0f s)' % (k, len(paths), time.time() - t0))
                sys.stdout.flush()
    return res


def gb(paths):
    return sum(os.path.getsize(p) for p in paths) / 1e9


def make(folders, manifest, workers, force):
    base = os.path.dirname(os.path.abspath(manifest))
    for f in folders:
        if not os.path.isdir(f):
            sys.exit('ERROR: folder not found: %s' % f)
        if same_or_inside(os.path.dirname(os.path.abspath(manifest)), f):
            sys.exit('ERROR: put the manifest outside the data folders (e.g. next to them)')
    if os.path.exists(manifest) and not force:
        sys.exit(
            'ERROR: %s already exists. It is your reference fingerprint, so it is never overwritten '
            '(use --force only if you really want a new reference).' % manifest
        )
    paths = []
    for f in folders:
        p = list_files(f)
        if not p:
            sys.exit('ERROR: no files found in %s - check the folder name. Nothing was written.' % f)
        ext = {}
        for q in p:
            e = os.path.splitext(q)[1] or '(none)'
            ext[e] = ext.get(e, 0) + 1
        print(
            '%-28s %7d files  %6.2f GB   %s'
            % (f, len(p), gb(p), ', '.join('%s %d' % kv for kv in sorted(ext.items())))
        )
        paths += p
    print('hashing %d files ...' % len(paths))
    h = hash_all(paths, workers)
    tmp = manifest + '.part'
    with open(tmp, 'w') as fh:
        for p in paths:
            fh.write('%s  %s\n' % (h[p], os.path.relpath(os.path.abspath(p), base).replace(os.sep, '/')))
    os.replace(tmp, manifest)
    with open(manifest + '.folders', 'w') as fh:
        for f in folders:
            fh.write(os.path.relpath(os.path.abspath(f), base).replace(os.sep, '/') + '\n')
    print('SHA-256 manifest written: %s (%d files)' % (manifest, len(paths)))


def verify(manifest, workers):
    if not os.path.isfile(manifest):
        sys.exit('ERROR: manifest not found: %s' % manifest)
    base = os.path.dirname(os.path.abspath(manifest))
    ref = {}
    with open(manifest) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line:
                h, rel = line.split('  ', 1)
                ref[os.path.normpath(os.path.join(base, rel))] = h
    if not ref:
        sys.exit('ERROR: the manifest lists no files')
    # folders that held the originals (saved by 'make'); used to find files added since then
    if os.path.isfile(manifest + '.folders'):
        with open(manifest + '.folders') as fh:
            roots = sorted(os.path.normpath(os.path.join(base, l.strip())) for l in fh if l.strip())
    else:
        roots = sorted(set(os.path.dirname(os.path.dirname(p)) for p in ref))
    present = [p for p in ref if os.path.isfile(p)]
    missing = sorted(set(ref) - set(present))
    new = sorted(
        set(os.path.normpath(q) for r in roots if os.path.isdir(r) for q in list_files(r)) - set(ref)
    )
    print('checking %d files in %s ...' % (len(ref), ', '.join(os.path.relpath(r) for r in roots)))
    h = hash_all(present, workers)
    changed = sorted(p for p in present if h[p] != ref[p])
    same = len(present) - len(changed)
    print(
        '\nchecked %d files: %d identical, %d changed, %d missing, %d new'
        % (len(ref), same, len(changed), len(missing), len(new))
    )
    for label, lst in (('CHANGED', changed), ('MISSING', missing), ('NEW', new)):
        for p in lst[:20]:
            print('   %-8s %s' % (label, os.path.relpath(p)))
        if len(lst) > 20:
            print('   ... and %d more %s files' % (len(lst) - 20, label.lower()))
    if changed or missing or new:
        print('RESULT: DIFFERENCES FOUND - see the list above')
        return 1
    print('RESULT: the original data are unchanged (all %d files byte-for-byte identical)' % same)
    return 0


def add_arguments(ap):
    ap.add_argument('action', choices=['make', 'verify'])
    ap.add_argument('folders', nargs='*', help='original data folders (make only)')
    ap.add_argument('-m', '--manifest', required=True, help='checksum file, e.g. originals_sha256.txt')
    ap.add_argument('--workers', type=int, default=min(8, os.cpu_count() or 2), help='parallel readers')
    ap.add_argument('--force', action='store_true', help='make: overwrite an existing manifest')


def run(a):
    if a.action == 'make':
        if not a.folders:
            sys.exit('ERROR: make needs the data folders, e.g. data/raw_sac data/raw_mseed')
        make(a.folders, a.manifest, a.workers, a.force)
    else:
        sys.exit(verify(a.manifest, a.workers))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(ap)
    run(ap.parse_args(argv))


if __name__ == '__main__':
    main()
