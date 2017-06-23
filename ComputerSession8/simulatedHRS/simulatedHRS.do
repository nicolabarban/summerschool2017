clear all
set maxvar 10000
use rand_matched, clear

keep subjectid raracem ragender rabyear raedyrs r1bmi  r1iwendy  r1height


gen age1=r1iwendy - rabyear

merge 1:1 subjectid  using profile_height
exit