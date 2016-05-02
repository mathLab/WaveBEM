At the top of the script.sh, there is a line containing all the FR
numbers we want to test for, and a base.prm file which is used 
for all the tests. 
The script creates the correct waveBem.prm in the correct directory, 
by substituting XXXX with the FR number. That FR number should be FR*1000, so 
.25 is 0250 (padded with zeros, to maintain alignment).

The script creates also a file "launch.sh", which you can run by

sh launch.sh

from hg1.

Step by step:

- edit base.prm
- edit script.sh (only the list of FR)
- launch script.sh (./script.sh)
- launch launch.sh (sh launch.sh)
