from lib.field import Field
from lib.util import *

m = '10011'
mod_length = len(m)
field = Field(int(m, 2))
elements = field.get_all_elems()
bin_elements = list(map(lambda it: tobin(it).zfill(mod_length - 1), elements))
for i, item in enumerate(bin_elements):
    print("%s: %s" % (str(i).zfill(2), item))