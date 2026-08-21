"""Expose the stable-boundary-layer mixing constants of VDFEXCU as namelist NAMVDFS.

WHY.  djf_bias_vertical_structure.py showed the coupled DJF land cold bias is
surface-confined and sits ENTIRELY above the screen: at 60-90N the low-level inversion
bias (T925-T2m) is +2.81 K while the surface inversion bias (T2m-Tskt) is +0.11 K.  The
skin is right; the 2m-925 hPa layer is not.  That layer is mixed by the stable branch of
VDFEXCU, whose controlling constants were all hardcoded literals, which is why every
surface-side lever tried so far (ECE_LAMSK_SN, Raupach z0, F1 RVZ0H) scored null -- they
adjust a part of the column that carries no error.

WHAT.  Five constants move from literals in VDFEXCU to module variables in YOEVDFS, set
in SUVDFS from the new namelist NAMVDFS:

    RSBLB     5.0    coefficient b of the LTG stable stability functions.  LOWER = longer
                     tail = more heat mixed down = weaker inversion = warmer winter screen
    RSBLD     1.0    coefficient d of the same functions
    RSBLLMIN  30.0   asymptotic mixing length above the stable BL, and floor within it [m]
    RSBLLMAX  300.0  cap on the stable mixing length [m]
    RSBLPBLF  0.1    fraction of the stable BL depth taken as the mixing length

Defaults equal the removed literals, so a fort.4 without &NAMVDFS is bit-identical to the
as-released model.  ZCB/ZCD/ZLMIN are used ONLY inside the IF (ZRI > 0) stable branch, so
these knobs cannot reach unstable or convective mixing -- checked before patching.

The namelist is read with a plain READ + IOSTAT rather than POSNAM, because POSNAM aborts
when the group is absent and no existing fort.4 carries &NAMVDFS.  IOSTAT < 0 (group not
present) leaves the defaults; IOSTAT > 0 (malformed group) ABORTS rather than being
swallowed -- silently ignoring a mistyped namelist is how this campaign has repeatedly
produced arms identical to their control.  SUVDFS echoes the values to NULOUT
unconditionally, so NODE.001_01 proves what the run actually used.

USAGE:  python expose_sbl_mixing.py <path-to-ifs-source>
Idempotent: re-running on an already-patched tree exits without changing anything.
"""
import sys
import os

SRC = sys.argv[1] if len(sys.argv) > 1 else sys.exit(__doc__)


def patch(rel, edits, marker):
    p = os.path.join(SRC, rel)
    s = open(p).read()
    if marker in s:
        print(f'{rel}: already patched, skipped')
        return
    for old, new in edits:
        if s.count(old) != 1:
            sys.exit(f'{rel}: expected exactly one match for:\n{old}\n'
                     f'found {s.count(old)} -- refusing to patch')
        s = s.replace(old, new)
    open(p, 'w').write(s)
    print(f'{rel}: patched')


patch('arpifs/module/yoevdfs.F90', [
    ("REAL(KIND=JPRB) :: DRI26",
     """REAL(KIND=JPRB) :: DRI26

!     Stable-boundary-layer mixing, exposed for tuning (namelist NAMVDFS).
!     These were literals inside VDFEXCU; the defaults set in SUVDFS reproduce
!     them exactly, so an unset namelist is bit-identical to the as-released code.
REAL(KIND=JPRB) :: RSBLB
REAL(KIND=JPRB) :: RSBLD
REAL(KIND=JPRB) :: RSBLLMIN
REAL(KIND=JPRB) :: RSBLLMAX
REAL(KIND=JPRB) :: RSBLPBLF"""),
    ("!     *RCHBA*     REAL       *CONSTANT A IN *HOLTSLAG AND *DEBRUIN",
     """!     *RSBLB*     REAL       *COEFFICIENT b OF THE LTG STABILITY FUNCTIONS IN
!                            STABLE CONDITIONS (VDFEXCU).  LOWERING IT LENGTHENS
!                            THE TAIL, MIXES MORE HEAT DOWN AND WEAKENS THE
!                            LOW-LEVEL INVERSION.  DEFAULT 5.
!     *RSBLD*     REAL       *COEFFICIENT d OF THE SAME FUNCTIONS.  DEFAULT 1.
!     *RSBLLMIN*  REAL       *ASYMPTOTIC MIXING LENGTH ABOVE THE STABLE BL AND
!                            FLOOR WITHIN IT [m].  DEFAULT 30.
!     *RSBLLMAX*  REAL       *CAP ON THE STABLE MIXING LENGTH [m].  DEFAULT 300.
!     *RSBLPBLF*  REAL       *FRACTION OF THE STABLE BL DEPTH TAKEN AS THE
!                            ASYMPTOTIC MIXING LENGTH.  DEFAULT 0.1.
!     *RCHBA*     REAL       *CONSTANT A IN *HOLTSLAG AND *DEBRUIN"""),
], marker='RSBLPBLF')

patch('arpifs/phys_ec/suvdfs.F90', [
    (""" & RCDHALF  ,RCDHPI2  ,RIMAX    ,DRITBL   ,DRI26  """,
     """ & RCDHALF  ,RCDHPI2  ,RIMAX    ,DRITBL   ,DRI26    ,&
 & RSBLB    ,RSBLD    ,RSBLLMIN ,RSBLLMAX ,RSBLPBLF
USE YOMLUN   , ONLY : NULOUT   ,NULNAM"""),
    ("INTEGER(KIND=JPIM) :: ITMAX, JIT, JJP",
     "INTEGER(KIND=JPIM) :: ITMAX, JIT, JJP\nINTEGER(KIND=JPIM) :: IVDFSIOS"),
    # NAMELIST is a specification statement and must precede the statement functions.
    ('#include "fcvdfs.func.h"',
     """NAMELIST /NAMVDFS/ RSBLB, RSBLD, RSBLLMIN, RSBLLMAX, RSBLPBLF

#include "abor1.intfb.h"
#include "fcvdfs.func.h\""""),
    ("""RCDHALF = 16._JPRB
RCDHPI2 = 2.0_JPRB*ATAN(1.0_JPRB)""",
     """RCDHALF = 16._JPRB
RCDHPI2 = 2.0_JPRB*ATAN(1.0_JPRB)

!     1.4 STABLE-BL MIXING IN VDFEXCU (namelist NAMVDFS)
!     Defaults are the values that were hardcoded in VDFEXCU, so leaving NAMVDFS
!     out of fort.4 reproduces the as-released model bit for bit.

RSBLB    = 5.0_JPRB
RSBLD    = 1.0_JPRB
RSBLLMIN = 30.0_JPRB
RSBLLMAX = 300.0_JPRB
RSBLPBLF = 0.1_JPRB

!     POSNAM is not used here: it ABORTS when the group is absent, and every
!     existing fort.4 lacks &NAMVDFS.  A plain namelist read distinguishes the two
!     cases that matter -- IOSTAT < 0 is "group not present", which must leave the
!     defaults standing, while IOSTAT > 0 is a malformed group, which must NOT be
!     swallowed.  Silently ignoring a mistyped namelist is how this configuration
!     has previously produced runs that were identical to their control.

REWIND(NULNAM)
READ(NULNAM,NML=NAMVDFS,IOSTAT=IVDFSIOS)
REWIND(NULNAM)
IF (IVDFSIOS > 0) THEN
  WRITE(NULOUT,'(A,I0)') ' SUVDFS: malformed &NAMVDFS, IOSTAT=',IVDFSIOS
  CALL ABOR1('SUVDFS: MALFORMED NAMVDFS')
ENDIF

!     Echoed unconditionally so the run log carries the values actually used.
WRITE(NULOUT,'(A,I0,5(1X,A,ES12.5))') ' SUVDFS: stable-BL mixing, iostat=', &
 & IVDFSIOS, &
 & 'RSBLB=',RSBLB,'RSBLD=',RSBLD,'RSBLLMIN=',RSBLLMIN, &
 & 'RSBLLMAX=',RSBLLMAX,'RSBLPBLF=',RSBLPBLF"""),
], marker='NAMVDFS')

patch('arpifs/phys_ec/vdfexcu.F90', [
    ("USE YOEVDFS   , ONLY : RCHBA, RCHBB, RCHBD, RCHB23A, RCHBBCD, RCHBCD, &",
     """USE YOEVDFS   , ONLY : RSBLB, RSBLD, RSBLLMIN, RSBLLMAX, RSBLPBLF
USE YOEVDFS   , ONLY : RCHBA, RCHBB, RCHBD, RCHB23A, RCHBBCD, RCHBCD, &"""),
    ("ZCD       = 1.0_JPRB", "ZCD       = RSBLD"),
    ("ZCB       = 5.0_JPRB", "ZCB       = RSBLB"),
    ("ZLMIN=30.0_JPRB", "ZLMIN=RSBLLMIN"),
    ("""          ZKLENT(JL,JK)=MAX(ZLMIN,ZP(JL)*0.1_JPRB*ZPBLHEIGHT(JL)*RPLRG)
        ELSE
          ZKLENT(JL,JK)=MAX(ZLMIN,0.1_JPRB*ZPBLHEIGHT(JL)*RPLRG)
        ENDIF
        ZKLENT(JL,JK)=MIN(300.0_JPRB,ZKLENT(JL,JK))""",
     """          ZKLENT(JL,JK)=MAX(ZLMIN,ZP(JL)*RSBLPBLF*ZPBLHEIGHT(JL)*RPLRG)
        ELSE
          ZKLENT(JL,JK)=MAX(ZLMIN,RSBLPBLF*ZPBLHEIGHT(JL)*RPLRG)
        ENDIF
        ZKLENT(JL,JK)=MIN(RSBLLMAX,ZKLENT(JL,JK))"""),
], marker='RSBLPBLF')
