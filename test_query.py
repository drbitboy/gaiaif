"""
Test query of Gaia star catalog SQLIte3 database (DB) files

Usage:

  python test_query.py

Setup:
  -  presence of gaia.sqlite3 and gaia_heavy.sqlite3 DB files or
     symlinks to same

"""

### Attach light and heavy DB files
attach1 = "attach 'gaia.sqlite3' as lgt"
attach2 = "attach 'gaia_heavy.sqlite3' as hvy"

### Select query with JOIN*:  ID offset*; RA; Dec; Source ID; Parallax.
### - WHERE clause limits returned rows to stars at high Declination
select1 = """
SELECT lgt.gaiartree.idoffset
     , lgt.gaialight.ra
     , lgt.gaialight.dec
     , hvy.gaiaheavy.source_id 
     , lgt.gaialight.parallax

FROM lgt.gaiartree

INNER JOIN lgt.gaialight
        ON lgt.gaialight.idoffset=lgt.gaiartree.idoffset

INNER JOIN hvy.gaiaheavy
        ON hvy.gaiaheavy.idoffset=lgt.gaiartree.idoffset

-- Use small range for Declination; this needs to be changed to 80
--   when the DB symlinks point to gaia_subset*.sqlite3 DB files
WHERE lgt.gaiartree.dechi > 89.95

ORDER BY hvy.gaiaheavy.source_id
;"""

import sqlite3 as sl3

### SQLite3 setup
cn = sl3.connect('')
cu = cn.cursor()
cu.execute(attach1)
cu.execute(attach2)

### Make the query
cu.execute(select1)

### Write results to STDOUT
print([ds[0] for ds in cu.description])
for row in cu: print(('{0:016x}'.format(row[-2]&0x7ffffff800000000),row,))
