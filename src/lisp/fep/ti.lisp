;;; Copyright (c) 2016, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :fep)

;;; What is this for????
(defparameter *vdw-bonded* "ifsc=1, scmask1=':1@H6', scmask2=':2@O1,H6'")

(defparameter *cando-charge-script*
  (let ((*package* (find-package :keyword)))
    (cl:format nil "~{~s~%~}"
               '((ql:quickload :fep)
                 (in-package :cando-user)
                 (defparameter *feps* (fep:load-feps ":%FEPS%"))
                 (read-am1-charges *feps*)
                 (calculate-am1-bcc-charges *feps*)
                 (fep:save-feps *feps* ":%OUTPUT%")
                 (core:exit)
                 ))))


(defparameter *solvate-addion-morph-side-script*
  (let ((*package* (find-package :keyword)))
    (cl:format nil "~{~s~%~}"
               '((ql:quickload :fep)
                 (in-package :cando-user)
                 (leap:setup-amber-paths)
                 (leap:source "leaprc.water.tip3p")
                 (leap:load-off "solvents.lib")
                 (leap:load-off "atomic_ions.lib")
                 (leap:load-atom-type-rules "ATOMTYPE_GFF.DEF")
                 (leap:source "leaprc.ff14SB.redq")
                 (leap:source "leaprc.gaff")
                 (defparameter *feps* (fep:load-feps ":%INPUT%"))
                 (defparameter *receptor* (first (fep:receptors *feps*)))
                 (defparameter *side-name* :%SIDE-NAME%)
                 (defparameter *morph* (find-morph-with-name :%MORPH-NAME% *feps*))
                 (defparameter *source* (fep:source *morph*))
                 (defparameter *target* (fep:target *morph*))
                 (fep:average-core-atom-positions *source* *target*)
                 (defparameter *ligands* (cando:combine (fep:molecule *source*)
                                                        (fep:molecule *target*)))
                 (format t "*side-name* --> ~a~%" *side-name*)
                 (cond
                   ((eq *side-name* :ligand)
                    (format t "Creating system by combining ligands~%")
                    (defparameter *system* (chem:matter-copy *ligands*)))
                   ((eq *side-name* :complex)
                    (format t "Creating system by combining ligands and receptor~%")
                    (defparameter *system* (cando:combine (chem:matter-copy *ligands*)
                                                          (chem:matter-copy *receptor*))))
                   (t (error "Unknown *side-name* ~s - must be one of :ligand or :complex" *side-name*)))
                 (leap:solvate-box *system*
                  (leap.core:lookup-variable (solvent-box *feps*))
                  (solvent-buffer *feps*)
                  (solvent-closeness *feps*))
                 (leap.add-ions:add-ions *system* :|Cl-| 0)
                 (cando:save-mol2 *system* (ensure-directories-exist ":%MOL2%"))
                 (ensure-directories-exist (pathname ":%TOPOLOGY%"))
                 (ensure-directories-exist (pathname ":%COORDINATES%"))
                 (leap.topology:save-amber-parm-format *system* ":%TOPOLOGY%" ":%COORDINATES%")
                 (core:exit)))))

(defparameter *prepare-min-in*
"minimisation
 &cntrl
   imin = 1, ntmin = 2,
   maxcyc = 100,
   ntpr = 20, ntwe = 20,
   ntb = 1,
   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, ifsc = 1, clambda = 0.5, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0,
   timask1 = ':%TIMASK1%', timask2 = ':%TIMASK2%',
   scmask1 = ':%SCMASK1%', scmask2 = ':%SCMASK2%'
 /
 &ewald
 / 
")

(defparameter *prepare-heat-in*
  "heating
 &cntrl
   imin = 0, nstlim = 10000, irest = 0, ntx = 1, dt = 0.002,
   nmropt = 1,
   ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
   tempi = 5.0, tautp = 1.0,
   vlimit = 20,
   ntb = 1,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, ifsc = 1, clambda = 0.5, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0,
   timask1 = ':%TIMASK1%', timask2 = ':%TIMASK2%',
   scmask1 = ':%SCMASK1%', scmask2 = ':%SCMASK2%'
 /
 &ewald
 / 

 &wt
   type='TEMP0',
   istep1 = 0, istep2 = 8000,                                      
   value1 = 5.0, value2 = 300.0
 /

 &wt type = 'END'
 /

")

(defparameter *prepare-press-in*
  "pressurising
 &cntrl
   imin = 0, nstlim = 10000, irest = 1, ntx = 5, dt = 0.002,
   ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
   tautp = 1.0,
   vlimit = 20,
   ntp = 1, pres0 = 1.0, taup = 2.0,
   ntb = 2,
   ntc = 2, ntf = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, ifsc = 1, clambda = 0.5, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0,
   timask1 = ':%TIMASK1%', timask2 = ':%TIMASK2%',
   scmask1 = ':%SCMASK1%', scmask2 = ':%SCMASK2%'
 /
 &ewald
 / 

")


(defparameter *cpptraj-strip*
"trajin :%COORDINATES%

# remove the two ligands and keep the rest
strip \":1,2\"
outtraj :%SOLVATED% onlyframes 1

# extract the first ligand
unstrip
strip \":2-999999\"
outtraj :%SOURCE% onlyframes 1

# extract the second ligand
unstrip
strip \":1,3-999999\"
outtraj :%TARGET% onlyframes 1
")

(defparameter *decharge*
  (let ((*package* (find-package :keyword)))
    (cl:format nil "~{~s~%~}"
               `((ql:quickload :fep)
                 (leap:setup-amber-paths)
                 (leap:source "leaprc.ff14SB.redq")
                 (leap:source "leaprc.gaff")
                 (leap:load-amber-params "frcmod.ionsjc_tip3p")
                 ;; load the fep-calculation
                 (defparameter *feps* (fep:load-feps ":%FEPS%"))
                 (defparameter lsolv (cando:load-mol2 ":%SOLVATED%"))
                 ;; after cpptraj molecules are merged and it drops solvent information - so we need to add it back
                 (cando:maybe-split-molecules-in-aggregate lsolv)
                 (cando:assign-solvent-molecules-using-residue-name lsolv (fep:solvent-residue-name *feps*))
                 (defparameter lsource (cando:load-mol2 ":%SOURCE%"))
                 ;;decharge transformation
                 (defparameter decharge (cando:combine (chem:matter-copy lsource)
                                                       (chem:matter-copy lsource)
                                                       (chem:matter-copy lsolv)))
                 (leap.set-box:set-box decharge :vdw)
                 (cando:save-mol2 decharge ":%DECHARGE-MOL2%")
                 (leap.topology:save-amber-parm-format decharge ":%DECHARGE-TOPOLOGY%" ":%DECHARGE-COORDINATES%")
                 (core:exit)
                 ))))

(defparameter *recharge*
  (let ((*package* (find-package :keyword)))
    (cl:format nil "~{~s~%~}"
               `((ql:quickload :fep)
                 (leap:setup-amber-paths)
                 (leap:source "leaprc.ff14SB.redq")
                 (leap:source "leaprc.gaff")
                 (leap:load-Amber-Params "frcmod.ionsjc_tip3p")
                 ;; load the fep-calculation
                 (defparameter *feps* (fep:load-feps ":%FEPS%"))
                 (defparameter lsolv (cando:load-mol2 ":%SOLVATED%"))
                 ;; after cpptraj molecules are merged and it drops solvent information - so we need to add it back
                 (cando:maybe-split-molecules-in-aggregate lsolv)
                 (cando:assign-solvent-molecules-using-residue-name lsolv (fep:solvent-residue-name *feps*))
                 (defparameter ltarget (cando:load-mol2 ":%TARGET%"))
                 (defparameter recharge (cando:combine (chem:matter-copy ltarget)
                                                       (chem:matter-copy ltarget)
                                                       (chem:matter-copy lsolv)))
                 (leap.set-box:set-box recharge :vdw)
                 (cando:save-mol2 recharge ":%RECHARGE-MOL2%")
                 (leap.topology:save-amber-parm-format recharge ":%RECHARGE-TOPOLOGY%" ":%RECHARGE-COORDINATES%")
                 (core:exit)
                 ))))

(defparameter *decharge-recharge-heat-in* 
  "heating
 &cntrl
   imin = 0, nstlim = 10000, irest = 0, ntx = 1, dt = 0.002,
   ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
   tempi = 50.0, tautp = 1.0,
   vlimit = 20,
   ntc = 2, ntf = 1,
   ntb = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

   nmropt = 1,
   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, clambda = :%LAMBDA%, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0,
   timask1 = ':%TIMASK1%', timask2 = ':%TIMASK2%',
   ifsc = :%IFSC%, CRGMASK = ':%CRGMASK%'
 /

 &ewald
 / 

 &wt
   type='TEMP0',
   istep1 = 0, istep2 = 8000,                                      
   value1 = 50.0, value2 = 300.0
 /

 &wt type = 'END'
 /
")

(defparameter *vdw-heat-in* 
  "heating
 &cntrl
   imin = 0, nstlim = 10000, irest = 0, ntx = 1, dt = 0.002,
   ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
   tempi = 50.0, tautp = 1.0,
   vlimit = 20,
   ntc = 2, ntf = 1,
   ntb = 1,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 1000, ntpr = 1000, ntwr = 5000,

   nmropt = 1,
   ntr = 1, restraint_wt = 5.00,
   restraintmask='!:WAT & !@H=',

   icfe = 1, clambda = :%LAMBDA%, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 0,
   timask1 = ':%TIMASK1%', timask2 = ':%TIMASK2%',
   ifsc = :%IFSC%, scmask1=':%SCMASK1%', scmask2=':%SCMASK2%', crgmask = ':%CRGMASK%'
 /

 &ewald
 / 

 &wt
   type='TEMP0',
   istep1 = 0, istep2 = 8000,
   value1 = 50.0, value2 = 300.0
 /

 &wt type = 'END'
 /
")

(defparameter *decharge-recharge-ti-in*
  "TI simulation
 &cntrl
   imin = 0, nstlim = 100000, irest = 1, ntx = 5, dt = 0.002,
   ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
   vlimit = 20,
   ntc = 2, ntf = 1,
   ntb = 2,
   ntp = 1, pres0 = 1.0, taup = 2.0,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 10000, ntpr = 10000, ntwr = 20000,

   icfe = 1, clambda = :%LAMBDA%, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 1,
   ifmbar = 1, bar_intervall = 1000, mbar_states = :%LAMBDA-WINDOWS-COUNT%,
   mbar_lambda = :%LAMBDA-WINDOWS%
   timask1 = ':%TIMASK1%', timask2 = ':%TIMASK2%',
   ifsc = :%IFSC%, crgmask = ':%CRGMASK%'
 /

 &ewald
 / 
")


(defparameter *vdw-ti-in*
  "TI simulation
 &cntrl
   imin = 0, nstlim = 100000, irest = 1, ntx = 5, dt = 0.002,
   ntt = 3, temp0 = 300.0, gamma_ln = 2.0, ig = -1,
   vlimit = 20,
   ntc = 2, ntf = 1,
   ntb = 2,
   ntp = 1, pres0 = 1.0, taup = 2.0,
   ioutfm = 1, iwrap = 1,
   ntwe = 1000, ntwx = 10000, ntpr = 10000, ntwr = 20000,

   icfe = 1, clambda = :%LAMBDA%, scalpha = 0.5, scbeta = 12.0,
   logdvdl = 1,
   ifmbar = 1, bar_intervall = 1000, mbar_states = :%LAMBDA-WINDOWS-COUNT%,
   mbar_lambda = :%LAMBDA-WINDOWS%
   timask1 = ':%TIMASK1%', timask2 = ':%TIMASK2%',
   ifsc = :%IFSC%, scmask1=':%SCMASK1%', scmask2=':%SCMASK2%', crgmask = ':%CRGMASK%'
 /

 &ewald
 / 
")


(defparameter *python-getdvdl*
  "import math
from collections import OrderedDict


class OnlineAvVar(object):
  '''A class that uses an online algorithm to compute mean and variance.'''

  def __init__(self, store_data = False):
    self.step = 0
    self.mean = 0.0
    self.M2 = 0.0

    self.store = store_data
    self.data = []


  def accumulate(self, x):
    '''Accumulate data points to compute mean and variance on-the-fly.'''

    self.step += 1

    delta = x - self.mean

    self.mean += delta / self.step
    self.M2 += delta * (x - self.mean)

    if self.store:
      self.data.append(x)


  def get_variance(self):
    '''Convenience funtion to return variance.'''

    return self.M2 / (self.step - 1)


  def get_stat(self):
    '''Convenience funtion to return mean and standard deviation.'''

    return self.mean, math.sqrt(self.M2 / (self.step - 1))

if __name__ == '__main__':
  import os, sys, glob
  import numpy as np

  prog = sys.argv[0]

  if len(sys.argv) < 4:
    print >>sys.stderr, 'Usage: %s skip glob_pattern windows' % prog
    sys.exit(1)

  skip = 5 # hard-coded - used to be passed on command line as int(sys.argv[1])
#  glob_pattern = sys.argv[2] # don't use glob_pattern
  output_filename = sys.argv[1]
  windows = sys.argv[2:]   # follows '--'
  extrap = 'polyfit' # or linear or polyfit
  stats = []

  data = OrderedDict()

  fout = open(output_filename,'w')
  for window in windows:
    cwd = os.getcwd()
#    os.chdir(os.path.dirname(os.path.realpath(window)))
###    print >>fout, 'window = %s' % window
    dVdl = OnlineAvVar()    
    ln = 0

    for en in [ window ]: # used to be glob.glob(glob_pattern):
###      print >>fout, 'en  %s' % en
      with open(en, 'r') as en_file:
        for line in en_file:
          ln += 1
          if ln > skip and line.startswith('L9') and not 'dV/dlambda' in line:
             val = line.split()[5]
###             print >>fout, 'ln = %d | val = %s' % (ln, val)
             dVdl.accumulate(float(val))

    mean, std = dVdl.get_stat()
    window_dir = os.path.dirname(os.path.realpath(window))
    window_name = window_dir[window_dir.rfind('/')+1:]
###    print >>fout, 'dVdl.step = %f' % dVdl.step
    data[float(window_name)] = (mean, std / math.sqrt(dVdl.step), std)

#    os.chdir(cwd)

  x = data.keys()
  y = [d[0] for d in data.values()]

  if extrap == 'linear':
    if 0.0 not in x:
      l = (x[0]*y[1] - x[1]*y[0]) / (x[0] - x[1])
      x.insert(0, 0.0)
      y.insert(0, l)

    if 1.0 not in x:
      l = ( (x[-2] - 1.0)*y[-1] + ((1.0-x[-1])*y[-2]) ) / (x[-2] - x[-1])
      x.append(1.0)
      y.append(l)
  elif extrap == 'polyfit' and (0.0 not in x or 1.0 not in x):
    if len(x) < 6:
      deg = len(x) - 1
    else:
      deg = 6

    coeffs = np.polyfit(x, y, deg)

    if 0.0 not in x:
      x.insert(0, 0.0)
      y.insert(0, coeffs[-1])

    if 1.0 not in x:
      x.append(1.0)
      y.append(sum(coeffs) )

  fout = open(output_filename,'w')
  print >>fout, '(( :lambdas ('
  for a, b in zip(x, y):
    if a in data:
      v = data[a]
      print >>fout, '( ', a, v[0], v[1], v[2], ') '
    else:
      print >> fout, ';;; ', a, b
  print >>fout, '))'
  print >>fout, '( :dg . ', np.trapz(y, x), '))'
  fout.close()

")


(defparameter *combine-stages*
  (let ((*package* (find-package :keyword)))
    (cl:format nil "~{~s~%~}"
               '((ql:quickload :fep)
                 (in-package :cando-user)
                 (let ((total 0.0)
		       parts)
                   (loop for filename in :%.parts%
                         for pathname = (pathname filename)
                         for contents = (with-open-file (fin pathname :direction :input) (read fin))
                         for dg = (cdr (assoc :dg contents))
                         do (format t "dg = ~f~%" dg)
                         (incf total dg)
			 (push (cons pathname dg) parts))
		   (let ((result (list (cons :%SIDE-NAME% (list (cons :parts parts) (cons :total total))))))
                     (with-open-file (fout ":%SIDE-ANALYSIS%" :direction :output)
				     (format fout "~s~%" result))))
                 (core:exit 0)))))

(defparameter *combine-sides*
  (let ((*package* (find-package :keyword)))
    (cl:format nil "~{~s~%~}"
               '((ql:quickload :fep)
                 (in-package :cando-user)
		 (let* ((parts
			 (loop for filename in :%.parts%
                               for pathname = (pathname filename)
                               for contents = (with-open-file (fin pathname :direction :input) (read fin))
                               append contents))
			(complex-parts (cdr (assoc :complex parts)))
			(_ (format t "complex-parts: ~s~%" complex-parts))
			(complex-total (cdr (assoc :total complex-parts)))
			(ligand-parts (cdr (assoc :ligand parts)))
			(_ (format t "ligand-parts: ~s~%" ligand-parts))
			(ligand-total (cdr (assoc :total ligand-parts)))
			(ddg (- complex-total ligand-total)) ; is this the right way around?
			(result (list (cons :%MORPH-NAME% (list* (cons :morph-delta-delta-g ddg)
								 parts)))))
                   (with-open-file (fout ":%MORPH-ANALYSIS%" :direction :output)
				   (format fout ";;; Is the sign right? complex-total - ligand-total~%")
				   (format fout "~s~%" result))
		   (format t "~s~%" result))
                 (core:exit 0)))))

(defclass ti-calculation ()
  ((compounds :initarg :compounds :accessor compounds)
   (paths :initarg :paths :accessor paths)))

(defclass ti-compound ()
  ((name :initarg :name :accessor name)))

(defun make-ti-compound (name)
  (make-instance 'ti-compound :name name))

;;; ------------------------------------------------------------
;;;
;;; File node types
;;;
;;; side means either :ligand or :complex
;;; lambda means the lambda value from 0.0 to 1.0
;;; morph refers to the morph between two compounds of a compound graph

(defclass node-file ()
  ((definers :initarg :definers :initform nil :accessor definers)
   (users :initarg :users :initform nil :accessor users)
   (name :initarg :name :accessor name)
   (extension :initarg :extension :accessor extension)))

(defclass feps-mixin () ())

(defclass feps-file (node-file feps-mixin)
  ()
  (:default-initargs
   :name "feps"
   :extension "cando"))

(defclass python-script-file (node-file)
  ((script :initarg :script :accessor script))
  (:default-initargs
   :name "getdvdl"
   :extension "py"))

(defclass cando-script-file (node-file)
  ((script :initarg :script :accessor script))
  (:default-initargs
   :name "cando"
   :extension "lisp"))

(defclass sqm-file-mixin ()
  ())

(defclass sqm-file (node-file sqm-file-mixin)
  ())

(defclass sqm-input-file (sqm-file)
  ()
  (:default-initargs
   :extension "in"))

(defclass sqm-atom-order-file (sqm-file)
  ()
  (:default-initargs
   :extension "order"))

(defclass sqm-output-file (sqm-file)
  ()
  (:default-initargs
   :extension "out"))

(defclass morph-file (node-file)
  ((morph :initarg :morph :accessor morph))
  (:documentation "Files that only depend on morph"))

(defclass morph-script (morph-file)
  ((script :initarg :script :accessor script)))

(defclass morph-side-file (morph-file)
  ((side :initarg :side :accessor side))
  (:documentation "Files that only depend on morph and side"))


(defclass morph-side-unknown-file (morph-side-file)
  ()
  (:documentation "Files with no specific purpose that depend on the morph/side"))

(defclass morph-side-mol2-file (morph-side-file)
  ()
  (:default-initargs
   :extension "mol2")
  (:documentation "mol2 files"))

(defclass morph-side-script (morph-side-file)
  ((script :initarg :script :accessor script))
  (:default-initargs
   :extension "in")
  (:documentation "This is the input script file"))

(defclass morph-side-restart-file (morph-side-file)
  ()
  (:default-initargs
   :extension "rst7")
  (:documentation "Restart files depending only on morph/side"))

(defclass morph-side-coordinate-file (morph-side-file)
  ()
  (:default-initargs
   :extension "crd")
  (:documentation "Coordinate files depending only on morph/side"))

(defclass morph-side-topology-file (morph-side-file)
  ()
  (:default-initargs
   :extension "parm7")
  (:documentation "Topology files depending only on morph/side"))

(defclass morph-side-trajectory-file (morph-side-file)
  ()
  (:documentation "Trajectory files depending only on morph/side"))

(defclass morph-side-stage-file (morph-side-file)
  ((stage :initarg :stage :accessor stage))
  (:documentation "Files that depend on morph, side, stage"))

(defclass morph-side-stage-script-file (morph-side-stage-file)
  ((script :initarg :script :accessor script))
  (:default-initargs
   :extension "in")
  (:documentation "Files that depend on morph, side, stage"))

(defclass morph-side-stage-mol2-file (morph-side-stage-file)
  ()
  (:default-initargs
   :extension "mol2"))

(defclass morph-side-stage-topology-file (morph-side-stage-file)
  ()
  (:default-initargs
   :extension "parm7"))

(defclass morph-side-stage-coordinates-file (morph-side-stage-file)
  ()
  (:default-initargs
   :extension "rst7"))

(defclass morph-side-stage-lambda-file (morph-side-stage-file)
  ((lambda% :initarg :lambda% :accessor lambda%))
  (:documentation "A type of file node identified by a side (ligand vs complex),
a lambda value and an morph." ))

(defclass morph-side-stage-lambda-unknown-file (morph-side-stage-lambda-file)
  ()
  (:documentation "A file that has an unknown purpose - use this as a placeholder until we know what
its for and then create a new class for it."))

(defclass morph-side-stage-lambda-amber-script (morph-side-stage-lambda-file)
  ((script :initarg :script :accessor script))
  (:default-initargs
   :extension "in")
  (:documentation "This is the input AMBER script file"))

(defclass morph-side-stage-lambda-trajectory-file (morph-side-stage-lambda-file)
  ()
  (:default-initargs
   :extension "crd")
  (:documentation "This is a restart file"))

(defclass morph-side-stage-lambda-restart-file (morph-side-stage-lambda-file)
  ()
  (:default-initargs
   :extension "rst7")
  (:documentation "This is a restart file"))

(defclass morph-side-stage-lambda-unknown-file (morph-side-stage-lambda-file)
  ()
  (:documentation "This is an unknown type file"))


(defgeneric node-pathname (node))

(defmethod node-pathname :around ((node node-file))
  (ensure-directories-exist (call-next-method)))

(defmethod node-pathname ((node feps-file))
  (make-pathname :name (string-downcase (name node))
                 :type (extension node)))

(defmethod node-pathname ((node node-file))
  (make-pathname :name (string-downcase (name node))
                 :type (extension node)))

(defmethod node-pathname ((node sqm-file))
  (merge-pathnames (call-next-method)
                   (make-pathname :directory (list
                                              :relative
                                              "am1bcc"))))

(defmethod node-pathname ((node morph-file))
  (merge-pathnames (make-pathname :name (string-downcase (name node)) :type (extension node))
                   (make-pathname :directory (list
                                              :relative
                                              (string (morph-string (morph node)))))))

(defmethod node-pathname ((node morph-side-file))
  (merge-pathnames (make-pathname :name (string-downcase (name node)) :type (extension node))
                   (make-pathname :directory (list
                                              :relative
                                              (string (morph-string (morph node)))
                                              (string-downcase (string (side node)))))))

(defmethod node-pathname ((node morph-side-stage-file))
  (merge-pathnames (make-pathname :name (string-downcase (name node)) :type (extension node))
                   (make-pathname :directory (list
                                              :relative
                                              (string (morph-string (morph node)))
                                              (string-downcase (string (side node)))
                                              (string-downcase (stage node))
                                              ))))

(defmethod node-pathname ((node morph-side-stage-lambda-file))
  (merge-pathnames (make-pathname :name (string-downcase (name node)) :type (extension node))
                   (make-pathname :directory (list
                                              :relative
                                              (string (morph-string (morph node)))
                                              (string-downcase (string (side node)))
                                              (string-downcase (stage node))
                                              (lambda% node)
                                              ))))

(defun option-extension (option)
  (case option
    (:-i "in")
    (:-o "out")
    (:-inf "info")
    (:-e "en")
    (:-r "rst7")
    (:-x "nc")
    (:-l "log")
    (otherwise (error "Unknown sander/pmemd option ~a" option))))


;;; ------------------------------------------------------------
;;;
;;;

(defclass job ()
  ((script :initform nil :initarg :script :accessor script)
   (inputs :initform nil :initarg :inputs :accessor inputs)
   (outputs :initform nil :initarg :outputs :accessor outputs)
   (parameters :initform nil :initarg :parameters :accessor parameters)
   (makefile-clause :initarg :makefile-clause :accessor makefile-clause)))

(defclass jupyter-job (job)
  ()
  (:default-initargs
   :makefile-clause nil))

(defclass cando-job (job)
  ())

(defclass read-charges-job (cando-job)
  ())
                     
(defclass scripted-job (job)
  ())

(defclass morph-job (scripted-job)
  ((morph :initarg :morph :accessor morph)))
  
(defclass sqm-job-mixin ()
  ())

(defclass cando-job-mixin ()
  ())

(defclass python-job-mixin ()
  ())

(defclass amber-job-mixin ()
  ())

(defclass cpptraj-job-mixin ()
  ())

(defclass morph-cando-job (morph-job cando-job-mixin)
  ())

(defclass ligand-sqm-job (job sqm-job-mixin)
  ((ligand-name :initarg :ligand-name :accessor ligand-name)))

(defclass morph-side-job (morph-job)
  ((side :initarg :side :accessor side))
  (:documentation "An AMBER job that only depends on the morph and side"))

(defclass morph-side-cando-job (morph-side-job cando-job-mixin)
  ())

(defclass morph-side-amber-job (morph-side-job amber-job-mixin)
  ())

(defclass morph-side-cpptraj-job (morph-side-job cpptraj-job-mixin)
  ())

(defclass morph-side-stage-job (morph-side-job)
  ((stage :initarg :stage :accessor stage))
  (:documentation "An AMBER job that only depends on the morph and side and stage"))

(defclass morph-side-stage-cando-job (morph-side-stage-job cando-job-mixin)
  ())

(defclass morph-side-stage-python-job (morph-side-stage-job python-job-mixin)
  ())

(defclass morph-side-stage-lambda-job (morph-side-stage-job)
  ((lambda% :initarg :lambda% :accessor lambda%)
   (lambda-values :initarg :lambda-values :accessor lambda-values))
  (:documentation "An AMBER job that only depends on the morph and side and stage and lambda"))

(defclass morph-side-stage-lambda-cando-job (morph-side-stage-lambda-job cando-job-mixin)
  ())

(defclass morph-side-stage-lambda-amber-job (morph-side-stage-lambda-job amber-job-mixin)
  ())



(defun job-substitutions (job)
  (let (substitutions)
    (push (cons (format nil "-i")
                (namestring (node-pathname (script job))))
          substitutions)
    (loop for input in (inputs job)
          for option = (option input)
          for file = (node input)
          if (eq option :.)
            do (let ((dot-parts (format nil "(list ~{\"~a\" ~})" (reverse (mapcar (lambda (x) (namestring (node-pathname x))) file)))))
                 (push (cons :%.PARTS% dot-parts) substitutions))
          else
            do (push (cons (intern (string-upcase option) :keyword)
                           (namestring (node-pathname file)))
                     substitutions))
    (loop for output in (outputs job)
          for option = (option output)
          for file = (node output)
          do (push (cons (intern (string-upcase option) :keyword)
                         (namestring (node-pathname file)))
                   substitutions))
    substitutions))

(defgeneric substitutions (calculation job node-file))

(defmethod substitutions (calculation job (node-file node-file))
  (job-substitutions job))

(defmethod substitutions (calculation job (node-file morph-file))
  (append
   (let* ((morph (morph node-file))
          (morph-mask (calculate-masks morph (mask-method calculation))))
     (list* (cons :%MORPH-NAME% (format nil "~s" (morph-string morph)))
            (mask-substitutions morph-mask)))
   (call-next-method)))

(defmethod substitutions (calculation job (node-file morph-side-file))
  (list* (cons "%SIDE-NAME%" (format nil "~s" (side node-file)))
         (call-next-method)))

(defmethod substitutions (calculation job (node-file morph-side-stage-file))
  (let* ((morph (morph node-file))
         (morph-mask (calculate-masks morph (mask-method calculation))))
    (append (list (cons "%STAGE-NAME%" (format nil "~s" (stage node-file))))
            (crgmask-substitutions morph-mask (stage node-file))
            (call-next-method))))

#+(or)
(defmethod substitutions (calculation (job t) (node-file morph-side-stage-lambda-file))
  (error "substitutions called with job ~s and node-file ~s" job node-file))

(defmethod substitutions (calculation (job morph-side-stage-lambda-amber-job) (node-file morph-side-stage-lambda-file))
  (append (list (cons :%LAMBDA% (format nil "~f" (lambda% node-file)))
                (cons :%LAMBDA-WINDOWS% (format nil "~{~f,~}" (lambda-values job)))
                (cons :%LAMBDA-WINDOWS-COUNT% (format nil "~a" (length (lambda-values job)))))
         (call-next-method)))

(defun output-file (amber-job option)
  (loop for output-arg in (outputs amber-job)
        when (eq option (option output-arg))
          do (return-from output-file (node output-arg)))
  (error "Could not find option ~a in ~a" option (outputs amber-job)))

(defclass step-info ()
  ((unique-step-name :initform (gensym) :accessor unique-step-name)))

(defun make-step-info (lam)
  (make-instance 'step-info :lam lam))

(defclass heat-step (step) ())

(defclass ti-step (step)
  ((lam :initarg :lam :accessor lam)))
(defclass ti-step-info (step-info)
  ((lam :initarg :lam :accessor lam)))

(defun make-ti-step-info (lam)
  (make-instance 'ti-step-info :lam lam))


(defun step-name (step)
  "This name is used to construct the directory for the step"
  (format nil "s-~,3f-~a" (lam step) (unique-step-name step)))

(defun make-morph-side-lambda-file (name extension)
  (make-ti-file :pathname (make-pathname :name (string name) :type (string extension))))

(defun connect-graph (amber-job)
  (let ((script (script amber-job))
        (inputs (inputs amber-job))
        (outputs (outputs amber-job)))
    (when script (push amber-job (users script)))
    (loop for input in inputs
          for node-input = (node input)
          if (consp node-input)
            do (mapc (lambda (one) (push amber-job (users one))) node-input)
          else
            do (push amber-job (users node-input)))
    (loop for output in outputs
          do (push amber-job (definers (node output))))
    amber-job))

(defclass argument ()
  ((option :initarg :option :accessor option)
   (node :initarg :node :accessor node)))

(defmethod print-object ((obj argument) stream)
  (print-unreadable-object (obj stream)
    (format stream "~a ~a" (option obj) (node obj))))

(defmethod users ((arg argument))
  (users (node arg)))

(defmethod definers ((arg argument))
  (definers (node arg)))

(defun arguments (&rest args)
  "Later we will add sanity checks to the args - for now just accumulate a list"
  (loop for cur = args then (cddr cur)
        for option = (first cur)
        for node = (second cur)
        until (null cur)
        while option
        collect (make-instance 'argument :option option :node node)))

(defun standard-makefile-clause (command)
  (format nil ":%OUTPUTS% : :%DEPENDENCY-INPUTS%
	runcmd -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	~a~%" command))

(defun standard-cando-makefile-clause (script &key add-inputs)
  (if add-inputs
      (standard-makefile-clause (format nil "icando-boehm -N -l ~a -- :%INPUTS%" (node-pathname script)))
      (standard-makefile-clause (format nil "icando-boehm -N -l ~a" (node-pathname script)))))

(defun make-ligand-sqm-step (ligand &key sqm-input-file)
  (connect-graph
   (make-instance 'ligand-sqm-job
                  :ligand-name (name ligand)
                  :inputs (arguments :-i sqm-input-file)
                  :outputs (arguments :-o (make-instance 'sqm-output-file :name (name ligand)))
                  :makefile-clause (standard-makefile-clause "sqm :%OPTION-INPUTS% -O :%OPTION-OUTPUTS%"))))

(defun make-morph-side-prepare-job (morph side &key executable name script input-coordinate-file input-topology-file)
  (unless executable (error "You must provide an executable"))
  (unless input-coordinate-file (error "You must provide an input-coordinate-file"))
  (unless input-topology-file (error "You must provide an input-topology-file"))
  (let ((script-file (make-instance 'morph-side-script :morph morph :side side :name name :script script)))
    (connect-graph
     (make-instance 'morph-side-amber-job
                    :morph morph
                    :side side
                    :script script-file
                    :inputs (arguments :-i script-file
                                       :-c input-coordinate-file ; min.rst7
                                       :-ref input-coordinate-file ; ../ligands_vdw_bonded.rst7
                                       :-p input-topology-file) ; $prmtop
                    :outputs (arguments :-o (make-instance 'morph-side-unknown-file
                                                           :morph morph :side side
                                                           :name name :extension "out")
                                        :-inf (make-instance 'morph-side-unknown-file
                                                             :morph morph :side side
                                                             :name name :extension "info")
                                        :-e (make-instance 'morph-side-unknown-file
                                                           :morph morph :side side
                                                           :name name :extension "en")
                                        :-r (make-instance 'morph-side-restart-file
                                                           :morph morph :side side
                                                           :name name)
                                        :-x (make-instance 'morph-side-trajectory-file
                                                           :morph morph :side side
                                                           :name name :extension "nc")
                                        :-l (make-instance 'morph-side-unknown-file
                                                           :morph morph :side side
                                                           :name name :extension "log"))
                    :makefile-clause (standard-makefile-clause (format nil "~a :%OPTION-INPUTS% -O :%OPTION-OUTPUTS%" executable))))))

(defun make-morph-side-prepare (morph side &rest args &key input-coordinate-file input-topology-file)
  (unless input-coordinate-file (error "You must provide an input-coordinate-file"))
  (unless input-topology-file (error "You must provide an input-topology-file"))
  (let* ((min-job (make-morph-side-prepare-job morph side
                                               :name "min"
                                               :executable "pmemd"
                                               :script *prepare-min-in*
                                               :input-coordinate-file input-coordinate-file
                                               :input-topology-file input-topology-file))
         (min.rst (output-file min-job :-r))
         (heat-job (make-morph-side-prepare-job morph side
                                                :name "heat"
                                                :executable "pmemd.cuda"
                                                :script *prepare-heat-in*
                                                :input-coordinate-file min.rst
                                                :input-topology-file input-topology-file))
         (heat.rst (output-file heat-job :-r))
         (press-job (make-morph-side-prepare-job morph side
                                                 :name "press"
                                                 :executable "pmemd.cuda"
                                                 :script *prepare-press-in*
                                                 :input-coordinate-file heat.rst
                                                 :input-topology-file input-topology-file)))
    press-job))


(defun make-morph-side-strip (morph side &key input-coordinate-file input-topology-file)
  (let ((script-file (make-instance 'morph-side-script :morph morph :side side :name "strip" :script *cpptraj-strip*)))
    (connect-graph
     (make-instance 'morph-side-cpptraj-job
                    :morph morph
                    :side side
                    :script script-file
                    :inputs (arguments :-i script-file :-p input-topology-file :%COORDINATES% input-coordinate-file)
                    :outputs (arguments :%SOLVATED% (make-instance 'morph-side-mol2-file
                                                                 :morph morph
                                                                 :side side
                                                                 :name "solvated"
                                                                 :extension "mol2")
                                        :%SOURCE% (make-instance 'morph-side-mol2-file
                                                               :morph morph
                                                               :side side
                                                               :name "source"
                                                               :extension "mol2")
                                        :%TARGET% (make-instance 'morph-side-mol2-file
                                                               :morph morph
                                                               :side side
                                                               :name "target"
                                                               :extension "mol2"))
                    :makefile-clause (standard-makefile-clause "cpptraj :%OPTION-INPUTS% :%OPTION-OUTPUTS%")))))

(defun make-heat-ti-step (morph side stage lam lambda-values &key input-coordinate-file input-topology-file)
  (let ((script (make-instance 'morph-side-stage-lambda-amber-script
                               :morph morph
                               :side side
                               :stage stage
                               :lambda% lam
                               :name "heat"
                               :script (cond
                                         ((eq stage :decharge)
                                          *decharge-recharge-heat-in*)
                                         ((eq stage :recharge)
                                          *decharge-recharge-heat-in*)
                                         ((eq stage :vdw-bonded)
                                          *vdw-heat-in*)
                                         ((null stage)
                                          *vdw-heat-in*)
                                         (t (error "Illegal stage ~a for make-heat-ti-step script - what do I do with this?" stage))))))
    (connect-graph
     (make-instance 'morph-side-stage-lambda-amber-job
                    :lambda% lam
                    :lambda-values lambda-values
                    :script script
                    :inputs (arguments :-i script
                                       :-c input-coordinate-file
                                       :-ref input-coordinate-file
                                       :-p input-topology-file)
                    :outputs (arguments :-o (make-instance 'morph-side-stage-lambda-unknown-file :morph morph :side side :stage stage :lambda% lam :name "heat" :extension "out")
                                        :-inf (make-instance 'morph-side-stage-lambda-unknown-file :morph morph :side side :stage stage :lambda% lam :name "heat" :extension "info")
                                        :-e (make-instance 'morph-side-stage-lambda-unknown-file :morph morph :side side :stage stage :lambda% lam :name "heat" :extension "en")
                                        :-r (make-instance 'morph-side-stage-lambda-restart-file :morph morph :side side :stage stage :lambda% lam :name "heat")
                                        :-x (make-instance 'morph-side-stage-lambda-trajectory-file :morph morph :side side :stage stage :lambda% lam :name "heat" :extension "nc")
                                        :-l (make-instance 'morph-side-stage-lambda-unknown-file :morph morph :side side :stage stage :lambda% lam :name "heat" :extension "log")
                                        )
                    :makefile-clause ":%OUTPUTS% : :%INPUTS%
	runcmd -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	pmemd.cuda -AllowSmallBox :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))))

(defun make-ti-step (morph side stage lam lambda-values &key input-coordinate-file input-topology-file)
  (let ((script (make-instance 'morph-side-stage-lambda-amber-script
                               :morph morph
                               :side side
                               :stage stage
                               :lambda% lam
                               :name "ti"
                               :script (cond
                                         ((eq stage :decharge) *decharge-recharge-ti-in*)
                                         ((eq stage :recharge) *decharge-recharge-ti-in*)
                                         ((eq stage :vdw-bonded) *vdw-ti-in*)
                                         ((null stage) *vdw-ti-in*)
                                         (t (error "Illegal stage ~a for make-ti-step script - what do I do with this?" stage))))))
    (connect-graph
     (make-instance 'morph-side-stage-lambda-amber-job
                    :lambda% lam
                    :lambda-values lambda-values
                    :script script
                    :inputs (arguments :-i script
                                       :-c input-coordinate-file
                                       :-p input-topology-file)
                    :outputs (arguments :-o (make-instance 'morph-side-stage-lambda-unknown-file :morph morph :side side :stage stage :lambda% lam :name "ti001" :extension "out")
                                        :-inf (make-instance 'morph-side-stage-lambda-unknown-file :morph morph :side side :stage stage :lambda% lam :name "ti001" :extension "info")
                                        :-e (make-instance 'morph-side-stage-lambda-unknown-file :morph morph :side side :stage stage :lambda% lam :name "ti001" :extension "en")
                                        :-r (make-instance 'morph-side-stage-lambda-restart-file :morph morph :side side :stage stage :lambda% lam :name "ti001")
                                        :-x (make-instance 'morph-side-stage-lambda-trajectory-file :morph morph :side side :stage stage :lambda% lam :name "ti001" :extension "nc")
                                        :-l (make-instance 'morph-side-stage-lambda-unknown-file :morph morph :side side :stage stage :lambda% lam :name "ti001" :extension "log")
                                        )
                    :makefile-clause ":%OUTPUTS% : :%INPUTS%
	runcmd -- :%DEPENDENCY-INPUTS% -- :%DEPENDENCY-OUTPUTS% -- \\
	pmemd.cuda -AllowSmallBox :%OPTION-INPUTS% \\
	  -O :%OPTION-OUTPUTS%"))))

(defclass ti-path ()
  ((lambdas :initform 11 :initarg :lambdas :accessor lambdas)
   (steps :initform nil :initarg :steps :accessor steps)
   (source-compound :initarg :source-compound :accessor source-compound)
   (target-compound :initarg :target-compound :accessor target-compound)))

(defparameter *identical-lambda-delta* 0.001)
(defmethod maybe-add-step (ti-path step)
  (prog1
      (pushnew step (steps ti-path) :test (lambda (x y)
                                            (and (eq (class-of x) (class-of y))
                                                 (let ((delta (abs (- (lam x) (lam y)))))
                                                   (< delta *identical-lambda-delta*)))))
    (format t "Number of steps ~a~%" (length (steps ti-path)))))

(defun replace-all (dict in-script)
  (let (script-result)
    (loop for script = in-script then script-result
          for (match . substitution) in dict
          for match-string = (cond
                               ((stringp match) (format nil ":~a" match))
                               ((symbolp match) (format nil "~s" match))
                               (t (error "Illegal match ~s for regex key" match)))
          do (setf script-result (cl-ppcre:regex-replace-all match-string script substitution))
          finally (return-from replace-all script-result))))

(defun makefile-substitutions (calculation job script-code)
  "This returns an alist of label/string substitutions.
The fancy part is the inputs - inputs that have the form :-xxx are added as option-inputs
If there is one (and only one) input with the argument :. - that is appended to the option-inputs with the
form '--' <list of file node pathnames>.  inputs and outputs with names like :yxxx where 'y' is NOT - or . are
added to inputs and outputs but not option-inputs or option-outputs"
  (let* (dependency-inputs
         inputs
         option-inputs
         outputs
         option-outputs
         (inputs-job (inputs job))
         (dot-option-arg (find-if (lambda (arg) (eq :. (option arg))) inputs-job)) ; Find argument with :. option
         (argument-inputs-job (remove-if (lambda (arg) (eq :. (option arg))) inputs-job)) ; Remove argument with :. option
         extra-substitutions)
    ;; Add the script as an input
    (when (script job)
      (push (namestring (node-pathname (script job))) dependency-inputs))
    ;; Add all inputs
    (loop for input in argument-inputs-job
          for node-input = (node input)
          do (pushnew (namestring (node-pathname node-input)) inputs :test #'string=)
             (if (char= (char (string (option input)) 0) #\-)
                 (progn
                   (push (string-downcase (option input)) option-inputs)
                   (push (namestring (node-pathname (node input))) option-inputs))
                 (push (cons (intern (string-upcase (option input)) :keyword) (namestring (node-pathname (node input)))) extra-substitutions)))
    (when dot-option-arg
      (mapc (lambda (one) (pushnew (namestring (node-pathname one)) inputs :test #'string=)) (reverse (node dot-option-arg)))
      (push "--" option-inputs)
      (mapc (lambda (one) (pushnew (namestring (node-pathname one)) option-inputs :test #'string=)) (reverse (node dot-option-arg))))
    (loop for output in (outputs job)
          do (pushnew (namestring (node-pathname (node output))) outputs :test #'string=)
             (if (char= (char (string (option output)) 0) #\-)
                 (progn
                   (push (string-downcase (option output)) option-outputs)
                   (push (namestring (node-pathname (node output))) option-outputs))
                 (push (cons (intern (string-upcase (option output)) :keyword) (namestring (node-pathname (node output)))) extra-substitutions)))
    (append
     (list* (cons :%DEPENDENCY-INPUTS% (format nil "~{~a ~}" (append dependency-inputs (reverse inputs))))
            (cons :%INPUTS% (format nil "~{~a ~}" (reverse inputs)))
            (cons :%DEPENDENCY-OUTPUTS% (format nil "~{~a ~}" (reverse outputs)))
            (cons :%OUTPUTS% (format nil "~{~a ~}" (reverse outputs)))
            (cons :%OPTION-INPUTS% (format nil "~{~a ~}" (reverse option-inputs)))
            (cons :%OPTION-OUTPUTS% (format nil "~{~a ~}" (reverse option-outputs)))
            (if script-code
                (list (cons :%SCRIPT-CODE% (format nil "~a" script-code)))
                nil))
     extra-substitutions)))

(defun write-file-if-it-has-changed (pathname code)
    (when (probe-file pathname)
      (let ((existing-code (with-open-file (fin pathname :direction :input)
                                    (let ((data (make-string (file-length fin))))
                                      (read-sequence data fin)
                                      data))))
        (when (string= code existing-code)
          (format t "Skipping generation of ~a - it has not changed~%" pathname)
          (return-from write-file-if-it-has-changed nil))))
    (format t "Generating script ~a~%" pathname)
    (with-open-file (fout (ensure-directories-exist pathname) :direction :output :if-exists :supersede)
      (write-string code fout)))

(defmethod generate-code (calculation (job job) makefile visited-nodes)
  ;; Generate script
  (let ((script (script job))
        script-code)
    (when script
      (let* ((raw-script (script script))
             (substituted-script (replace-all (substitutions calculation job script) raw-script)))
        (setf script-code substituted-script)
        (write-file-if-it-has-changed (node-pathname script) substituted-script)))
    (let ((raw-makefile-clause (makefile-clause job)))
      (when raw-makefile-clause
        (let* ((makefile-substitutions (makefile-substitutions calculation job script-code))
               (substituted-makefile-clause (replace-all makefile-substitutions raw-makefile-clause)))
          (write-string substituted-makefile-clause makefile))
        (terpri makefile)
        (terpri makefile)))
    (loop for output in (outputs job)
          do (loop for child in (users output)
                   unless (gethash child visited-nodes)
                     do (setf (gethash child visited-nodes) t)
                        (generate-code calculation child makefile visited-nodes)))))

(defun generate-all-code (calculation work-list final-outputs)
  (with-top-directory (calculation)
    (let ((visited-nodes (make-hash-table))
          (makefile-pathname (ensure-directories-exist (merge-pathnames "makefile"))))
      (format t "Writing makefile to ~a~%" (translate-logical-pathname makefile-pathname))
      (let ((body (with-output-to-string (makefile)
                    (loop for job in work-list
                          do (generate-code calculation job  makefile visited-nodes)))))
        (write-file-if-it-has-changed
         makefile-pathname
         (with-output-to-string (makefile)
           (format makefile "# cando version ~a~%" (lisp-implementation-version))
           (format makefile "# fep version ~a~%" (fep:version))
           (format makefile "all : ~{~a ~}~%" (mapcar (lambda (file) (node-pathname file)) final-outputs))
           (format makefile "~aecho DONE~%" #\tab)
           (format makefile "~%")
           (write-string body makefile)
           (terpri makefile)))))))

#|
(defparameter *morphs*
  (let* ((x1 (make-ti-compound "x1"))
         (x2 (make-ti-compound "x2"))
         (x3 (make-ti-compound "x3"))
         (path1 (make-ti-path 11 x1 x2))
         (_ (generate-steps path1))
         (path2 (make-ti-path 11 x2 x3))
         (_ (generate-steps path2))
         (calc (make-instance 'ti-calculation :compounds (list x1 x2 x3) :paths (list path1 path2))))
    calc))

*morphs*

(defparameter *morphs* (make-instance 'ti-calculation))
(cando:save-cando *morphs* "/tmp/calc")

(cando:save-cando (make-instance 'foo) "/tmp/test")


(defclass bar () ()) (let ((*print-readably* t)) (print (make-instance 'bar)))

(defmethod print-object ((x foo) stream)
  (format stream "#S(foo)"))

(generate-all-scripts *morphs*)

|#

(defmethod print-object ((object chem:aggregate) stream)
  "Aggregates can have atom graphs that are way too wide and deep to print the 
normal way - so we short circuit it here using a mol2 file"
  (let ((agg-string (chem:aggregate-as-mol2-string object)))
    (format stream "#.(with-input-from-string (sin ~s) (chem:read-mol2 sin)) " agg-string)))


(defun save-feps (feps feps-file)
  (cando:save-cando feps feps-file))

(defun load-feps (feps-file)
  (cando:load-cando feps-file))

#+(or)
(defun save-feps (feps feps-file)
  (setf (receptor-strings feps) nil)
  (let* ((receptors (receptors feps))
         reversed-receptor-strings)
    (loop for receptor in receptors
          for receptor-string = (chem:aggregate-as-mol2-string receptor)
          do (push receptor-string reversed-receptor-strings))
    (setf (receptor-strings feps) (nreverse reversed-receptor-strings)))
  (cando:save-cando feps feps-file))


#+(or)
(defun load-feps (feps-file)
  "Load a feps database and register the topologys"
  (let ((feps (cando:load-cando feps-file))
        reversed-receptors)
    (loop for receptor-string in (receptor-strings feps)
          for receptor = (with-input-from-string (sin receptor-string)
                           (chem:read-mol2 sin))
          do (push receptor reversed-receptors))
    (setf (receptors feps) (nreverse reversed-receptors))
    (cando:register-topology (chem:get-name (core-topology feps)) (core-topology feps))
    (maphash (lambda (part-name name-topologys)
               (format t "name-topologys = ~s~%" name-topologys)
               (loop for (name . topology) in name-topologys
                     do (format t "Registering ~a ~a~%" name topology)
                     do (cando:register-topology name topology)))
             (side-topologys feps))
    feps))

