attach1 = "attach 'gaia.sqlite3' as lgt"
attach2 = "attach 'gaia_heavy.sqlite3' as hvy"
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
WHERE lgt.gaiartree.dechi > 80.00

ORDER BY hvy.gaiaheavy.source_id
;"""

import sqlite3 as sl3

cn = sl3.connect('')
cu = cn.cursor()
cu.execute(attach1)
cu.execute(attach2)
cu.execute(select1)
print([ds[0] for ds in cu.description])
for row in cu: print(('{0:016x}'.format(row[-2]&0x7ffffff800000000),row,))
