(in-package :leap.add-ions)

(defun ion-topology-atom-type (ion-topology)
  (let* ((constitution (chem:get-constitution ion-topology))
         (constitution-atoms (chem:constitution/get-constitution-atoms constitution))
         (number-of-atoms (chem:number-of-atoms constitution-atoms)))
    (if (= number-of-atoms 1)
        (let ((constitution-atom (chem:atom-with-id constitution-atoms 0)))
          (chem:atom-type constitution-atom))
        (error "The topology ~s must have only one atom" ion-topology))))

(defun add-ions (mol ion1 ion1-number &optional ion2 ion2-number)
  (warn "This function is a bit broken - every molecule must specify what force-field it will use - leap.core::merged-force-field is incorrect - the force-field has to be looked up from the ion-topology since the ion topology will create single residue in a molecule.  Using :default for now")
  (let* ((force-field-name :default)
         (nonbond-db (leap.core:nonbond-force-field-component force-field-name))
         (ion1-topology (cando:lookup-topology ion1))
         (ion1-type (ion-topology-atom-type ion1-topology))
         (ion1-type-index (chem:find-type-index nonbond-db ion1-type))
         (ion1-ffnonbond (chem:get-ffnonbond-using-type-index nonbond-db ion1-type-index))
         (ion1-size (chem:get-radius-angstroms ion1-ffnonbond))
         (ion1-residue (chem:build-residue-single-name ion1-topology))
         (ion1-atom (chem:content-at ion1-residue 0))
         (ion1-mol (chem:make-molecule))
         (ion1-agg (chem:make-aggregate))
         (solvent-vec (make-array 100 :fill-pointer 0 :adjustable t))
         (include-solvent 0)
         (target-charge 0)
         (ion-min-size 0.0)
         (grid-space 1.0)
         (shell-extent 4.0)
         (at-octree 1)
         (dielectric 1)
         (ion2-size 0.0) 
         (ion2-type-index 0)
         ion2-ffnonbond)
    (chem:setf-force-field-name ion1-mol force-field-name)
    (chem:add-matter ion1-mol ion1-residue)
    (chem:add-matter ion1-agg ion1-mol)
    (if (and ion2 (= 0 ion2-number))
        (error "'0' is not allowed as the value for the second ion."))
    ;;Make array of solvent residues
    (chem:map-residues 
     nil
     (lambda (r)
       (when (eq (chem:get-name r) :wat)
         (chem:setf-residue-type r :solvent) 
         (vector-push-extend r solvent-vec)))
     mol)
    ;;Consider target unit's charge
    (chem:map-atoms
     nil
     (lambda (r)
       (setf target-charge (+ (chem:get-charge r) target-charge)))
     mol)
    (if (and (= 0 target-charge)
             (= 0 ion1-number))
        (progn
          (format t "The total charge is ~a and so there are no charges to neutralize.~%" target-charge)
          (return-from add-ions))
        (format t "Total charge ~,2f~%" target-charge))
    ;;Consider neutralizetion    
    (if (= ion1-number 0)
        (progn
          (when (or (and (< (chem:get-charge ion1-atom) 0)
                       (< target-charge 0))
                  (and (> (chem:get-charge ion1-atom) 0)
                       (> target-charge 0)))
            (warn "1st Ion and target aggregate have charges of the same sign - aggregate charge: ~f; ion charge: ~f"
                  target-charge (chem:get-charge ion1-atom))
            (return-from add-ions nil))
          ;;Get the nearest integer number of ions that we need to add to get as close as possible to neutral.
          (setf ion1-number (round (/ (abs target-charge) (abs (chem:get-charge ion1-atom)))))
          (when ion2
              (warn "Neutralization - can't do 2nd ion.")
              (return-from add-ions nil))
          (format t "~d ~a ions required to neutraize. ~%" ion1-number (chem:get-name ion1-atom))))
    ;;Consider ion sizes and postions
    (if ion2
        (let* ((ion2-topology (cando:lookup-topology ion2))
               (ion2-type (ion-topology-atom-type ion2-topology)))
          (setf ion2-type-index (chem:find-type-index nonbond-db ion2-type))
          (setf ion2-ffnonbond (chem:get-ffnonbond-using-type-index nonbond-db ion2-type-index))
          (setf ion2-size (chem:get-radius-angstroms ion2-ffnonbond))
          (setf ion-min-size (min ion1-size ion2-size)))
        (setf ion-min-size ion1-size))
    (format t "Adding ~d counter ions to ~a using 1A grid. ~%"
            (if ion2 (+ ion1-number ion2-number) ion1-number) (chem:get-name mol))
    (if ion2
        (if (> (+ ion1-number ion2-number) 5)
            (setf ion-min-size (if (> ion1-size ion2-size)
                                   (* ion1-size
                                      (if (> (exp (/ (log (+ ion1-number ion2-number 1.0)) 3.0)) 1.0)
                                          (exp (/ (log (+ ion1-number ion2-number 1.0)) 3.0))
                                          1.0))
                                   (* ion2-size
                                      (if (> (exp (/ (log (+ ion1-number ion2-number 1.0)) 3.0)) 1.0)
                                          (exp (/ (log (+ ion1-number ion2-number 1.0)) 3.0))
                                          1.0)))))
        (if (> ion1-number 5)
            (setf ion-min-size (* ion1-size
                                  (if (> (exp (/ (log (+ ion1-number 1.0)) 3.0)) 1.0)
                                      (exp (/ (log (+ ion1-number 1.0)) 3.0))
                                      1.0)))))
    
    ;;Build grid and calc potential on it
    (let ((octree (core:make-cxx-object 'chem:add-ion-octree))
          (ion2-mol (chem:make-molecule))
          (ion2-agg (chem:make-aggregate))
          ion2-topology ion2-residue ion2-atom)
      (if ion2
          (progn 
            (setf ion2-topology (cando:lookup-topology ion2)
                  ion2-residue (chem:build-residue-single-name ion2-topology)
                  ion2-atom (chem:content-at ion2-residue 0))
            (chem:add-matter ion2-mol ion2-residue)
            (chem:add-matter ion2-agg ion2-mol)))
      (format t "About to octree-create~%")
      (chem:octree-create octree mol :shell grid-space ion-min-size shell-extent nonbond-db include-solvent t)
      (format t "Came out of octree-create~%")
      (if (aref solvent-vec 0)
          (format t " Solvent present: replacing closest with ion when steric overlaps occur~%")
          (format t "(no solvent present)~%"))
      (multiple-value-bind (min-charge-point max-charge-point)
          (chem:octree-init-charges octree at-octree dielectric ion1-size)
        (loop 
           for new-point = (if (< (chem:get-charge ion1-atom) 0) max-charge-point min-charge-point)
           for ion1-copy = (chem:matter-copy ion1-agg)
           for ion1-transform = (geom:make-m4-translate new-point)
           while (if ion2
                     (or (> ion1-number 0)
                         (> ion2-number 0))
                     (> ion1-number 0))
           do (if (> ion1-number 0)
                  (progn
                    (if (aref solvent-vec 0)
                        (check-solvent mol solvent-vec ion1-copy new-point))  
                    (chem:apply-transform-to-atoms ion1-copy ion1-transform)
                    (chem:add-matter mol ion1-copy)
                    (format t "Placed ~a in ~a at (~,2f, ~,2f, ~,2f).~%" (chem:get-name ion1-atom) 
                            (chem:get-name mol)
                            (geom:vx new-point)
                            (geom:vy new-point)
                            (geom:vz new-point))
                    (chem:octree-delete-sphere octree new-point (if ion2
                                                                  (+ ion1-size ion2-size)
                                                                  (+ ion1-size ion1-size)))
                    (multiple-value-bind (min-new-charge-point max-new-charge-point)
                        (chem:octree-update-charge octree new-point (chem:get-charge ion1-atom)
                                                     (if ion2 ion2-size ion1-size))
                      (setf min-charge-point min-new-charge-point)
                      (setf max-charge-point max-new-charge-point))
                    (decf ion1-number)))
           do (if (and ion2 (> ion2-number 0))
                  (progn 
                    (if (< (chem:get-charge ion2-atom) 0)
                        (setf new-point max-charge-point)
                        (setf new-point min-charge-point))
                    (let ((ion2-copy (chem:matter-copy ion2-agg))
                          (ion2-transform (geom:make-m4-translate new-point)))
                      (chem:apply-transform-to-atoms ion2-copy ion2-transform)
                      (chem:add-matter mol ion2-copy)
                      (format t "Placed ~a in ~a at (~,2f, ~,2f, ~,2f).~%" (chem:get-name ion2-atom) 
                              (chem:get-name mol)
                              (geom:vx new-point)
                              (geom:vy new-point)
                              (geom:vz new-point))
                      (chem:octree-delete-sphere octree new-point (if ion1
                                                                          (+ ion2-size ion1-size)
                                                                          (+ ion2-size ion2-size)))
                      (multiple-value-bind (min-new-charge-point max-new-charge-point)
                          (chem:octree-update-charge octree new-point (chem:get-charge ion2-atom)
                                                       (if ion1 ion1-size ion2-size))
                        (setf min-charge-point min-new-charge-point)
                        (setf max-charge-point max-new-charge-point)))
                    (decf ion2-number))))))))


(defun add-ions-2 (mol ion1 ion1-number &optional ion2 ion2-number)
  (warn "This function is now broken - every molecule must specify what force-field it will use - leap.core::merged-force-field is incorrect - the force-field has to be looked up from the ion-topology since the ion topology will create single residue in a molecule. Using :default")
  (let* ((force-field-name :default)
         (nonbond-db (leap.core:nonbond-force-field-component force-field-name))
         (ion1-topology (cando:lookup-topology ion1))
         (ion1-type (ion-topology-atom-type ion1-topology))
         (ion1-type-index (chem:find-type-index nonbond-db ion1-type))
         (ion1-ffnonbond (chem:get-ffnonbond-using-type-index nonbond-db ion1-type-index))
         (ion1-size (chem:get-radius-angstroms ion1-ffnonbond))
         (ion1-residue (chem:build-residue-single-name ion1-topology))
         (ion1-atom (chem:content-at ion1-residue 0))
         (ion1-mol (chem:make-molecule))
         (ion1-agg (chem:make-aggregate))
         (include-solvent 1)
         (target-charge 0)
         (ion-min-size 0.0)
         (grid-space 1.0)
         (shell-extent 4.0)
         (at-octree 1)
         (dielectric 1)
         (ion2-size 0.0) 
         (ion2-type-index 0)
         ion2-ffnonbond)
    (chem:add-matter ion1-mol ion1-residue)
    (chem:add-matter ion1-agg ion1-mol)
    (if (and ion2 (= 0 ion2-number))
        (error "'0' is not allowed as the value for the second ion."))
    ;;Consider target unit's charge
    (chem:map-atoms
     nil
     (lambda (r)
       (setf target-charge (+ (chem:get-charge r) target-charge)))
     mol)
    (if (and (= 0 target-charge)
             (= 0 ion1-number))
        (progn
          (format t "The total charge is ~a and so there are no charges to neutralize.~%" target-charge)
          (return-from add-ions-2))
        (format t "Total charge ~f~%" target-charge))
    ;;Consider neutralizetion    
    (if (= ion1-number 0)
        (progn
          (when (or (and (< (chem:get-charge ion1-atom) 0)
                       (< target-charge 0))
                  (and (> (chem:get-charge ion1-atom) 0)
                       (> target-charge 0)))
            (warn "1st Ion and target aggregate have charges of the same sign - aggregate charge: ~f; ion charge: ~f"
                  target-charge (chem:get-charge ion1-atom))
            (return-from add-ions-2 nil))
          ;;Get the nearest integer number of ions that we need to add to get as close as possible to neutral.
          (setf ion1-number (round (/ (abs target-charge) (abs (chem:get-charge ion1-atom)))))
          (when ion2
              (error "Neutralization - can't do 2nd ion."))
          (format t "~d ~a ions required to neutraize. ~%" ion1-number (chem:get-name ion1-atom))))
    ;;Consider ion sizes and postions
    (if ion2
        (let* ((ion2-topology (cando:lookup-topology ion2))
               (ion2-type (ion-topology-atom-type ion2-topology)))
          (setf ion2-type-index (chem:find-type-index nonbond-db ion2-type))
          (setf ion2-ffnonbond (chem:get-ffnonbond-using-type-index nonbond-db ion2-type-index))
          (setf ion2-size (chem:get-radius-angstroms ion2-ffnonbond))
          (setf ion-min-size (min ion1-size ion2-size)))
        (setf ion-min-size ion1-size))
    (format t "Adding ~d counter ions to ~a using 1A grid. ~%"
            (if ion2 (+ ion1-number ion2-number) ion1-number) (chem:get-name mol))
    (if ion2
        (if (> (+ ion1-number ion2-number) 5)
            (setf ion-min-size (if (> ion1-size ion2-size)
                                   (* ion1-size
                                      (if (> (exp (/ (log (+ ion1-number ion2-number 1.0)) 3.0)) 1.0)
                                          (exp (/ (log (+ ion1-number ion2-number 1.0)) 3.0))
                                          1.0))
                                   (* ion2-size
                                      (if (> (exp (/ (log (+ ion1-number ion2-number 1.0)) 3.0)) 1.0)
                                          (exp (/ (log (+ ion1-number ion2-number 1.0)) 3.0))
                                          1.0)))))
        (if (> ion1-number 5)
            (setf ion-min-size (* ion1-size
                                  (if (> (exp (/ (log (+ ion1-number 1.0)) 3.0)) 1.0)
                                      (exp (/ (log (+ ion1-number 1.0)) 3.0))
                                      1.0)))))
    
    ;;Build grid and calc potential on it
    (let ((octree (core:make-cxx-object 'chem:add-ion-octree))
          (ion2-mol (chem:make-molecule))
          (ion2-agg (chem:make-aggregate))
          ion2-topology ion2-residue ion2-atom)
      (if ion2
          (progn 
            (setf ion2-topology (cando:lookup-topology ion2)
                  ion2-residue (chem:build-residue-single-name ion2-topology)
                  ion2-atom (chem:content-at ion2-residue 0))
            (chem:add-matter ion2-mol ion2-residue)
            (chem:add-matter ion2-agg ion2-mol)))
      (format t "About to octree-create~%")
      (chem:octree-create octree mol :shell grid-space ion-min-size shell-extent nonbond-db include-solvent t)
      (format t "Came out of octree-create~%")
      (multiple-value-bind (min-charge-point max-charge-point)
          (chem:octree-init-charges octree at-octree dielectric ion1-size)
        (loop 
           for new-point = (if (< (chem:get-charge ion1-atom) 0) max-charge-point min-charge-point)
           for ion1-copy = (chem:matter-copy ion1-agg)
           for ion1-transform = (geom:make-m4-translate new-point)
           while (if ion2
                     (or (> ion1-number 0)
                            (> ion2-number 0))
                     (> ion1-number 0))
           do (if (> ion1-number 0)
                  (progn
                    (chem:apply-transform-to-atoms ion1-copy ion1-transform)
                    (chem:add-matter mol ion1-copy)
                    (format t "Placed ~a in ~a at (~,2f, ~,2f, ~,2f).~%" (chem:get-name ion1-atom) 
                            (chem:get-name mol)
                            (geom:vx new-point)
                            (geom:vy new-point)
                            (geom:vz new-point))
                    (chem:octree-delete-sphere octree new-point (if ion2
                                                                        (+ ion1-size ion2-size)
                                                                        (+ ion1-size ion1-size)))
                    (multiple-value-bind (min-new-charge-point max-new-charge-point)
                        (chem:octree-update-charge octree new-point (chem:get-charge ion1-atom)
                                                     (if ion2 ion2-size ion1-size))
                      (setf min-charge-point min-new-charge-point)
                      (setf max-charge-point max-new-charge-point))
                    (decf ion1-number)))
           do (if (and ion2 (> ion2-number 0))
                  (progn 
                    (if (< (chem:get-charge ion2-atom) 0)
                        (setf new-point max-charge-point)
                        (setf new-point min-charge-point))
                    (let ((ion2-copy (chem:matter-copy ion2-agg))
                          (ion2-transform (geom:make-m4-translate new-point)))
                      (chem:apply-transform-to-atoms ion2-copy ion2-transform)
                      (chem:add-matter mol ion2-copy)
                      (format t "Placed ~a in ~a at (~,2f, ~,2f, ~,2f).~%" (chem:get-name ion2-atom) 
                              (chem:get-name mol)
                              (geom:vx new-point)
                              (geom:vy new-point)
                              (geom:vz new-point))
                      (chem:octree-delete-sphere octree new-point (if ion1
                                                                          (+ ion2-size ion1-size)
                                                                          (+ ion2-size ion2-size)))
                      (multiple-value-bind (min-new-charge-point max-new-charge-point)
                          (chem:octree-update-charge octree new-point (chem:get-charge ion2-atom)
                                                       (if ion1 ion1-size ion2-size))
                        (setf min-charge-point min-new-charge-point)
                        (setf max-charge-point max-new-charge-point)))
                    (decf ion2-number))))))))


(defun check-solvent (mol solvent-vec ion1-copy new-point)
  (let ((dmin2 100000000)
        (closest-atom-vector (make-array 3 :element-type 'double-float))
        (d2 0)
        (x 0)
        (y 0)
        (z 0))
    (loop for i from 0 below (length solvent-vec)
          for residue = (aref solvent-vec i)
          do (chem:map-atoms
              nil
              (lambda (a)
                (setf x (- (geom:vx new-point) (geom:vx (chem:get-position a))))
                (setf y (- (geom:vy new-point) (geom:vy (chem:get-position a))))
                (setf z (- (geom:vz new-point) (geom:vz (chem:get-position a))))
                (setf d2 (+ (* x x) (* y y) (* z z)))
                (if (< d2 dmin2)
                    (progn
                      (setf dmin2 d2)
                      (setf (aref closest-atom-vector 0) (geom:vx (chem:get-position a)))
                      (setf (aref closest-atom-vector 1) (geom:vy (chem:get-position a)))
                      (setf (aref closest-atom-vector 2) (geom:vz (chem:get-position a))))))
              residue))
     (if (< dmin2 9)
        (progn 
          (chem:map-residues 
           nil
           (lambda (r)
             (chem:map-atoms
              nil
              (lambda (a)
                (when (and (= (geom:vx (chem:get-position a)) (aref closest-atom-vector 0))
                           (= (geom:vy (chem:get-position a)) (aref closest-atom-vector 1))
                           (= (geom:vz (chem:get-position a)) (aref closest-atom-vector 2)))
                  (chem:set-name r :delete)))
              r))
           mol)
          ;; If molecule has one water residue, delete molecule.
          ;; If molecule has more than two residues, delete one residue.
          (chem:map-molecules
           nil
           (lambda (m)
             (if (= (chem:residue-count m) 1)
                   (when (eq (chem:get-name (chem:get-residue m 0)) :delete)
                     (format t "Replacing solvent molecule~%")
                     (chem:remove-molecule mol m))
                   (chem:map-residues
                    nil
                    (lambda (r)
                      (when (eq (chem:get-name r) :delete)
                        (format t "Replacing solvent molecule~%")
                        (chem:remove-residue m r)))
                    m)))
           mol)))))

(defun add-ions-rand (aggregate ion1 ion1-number &key ion2 ion2-number (separation 0.0))
  (let* ((energy-function (chem:make-energy-function :matter aggregate
                                                     :use-excluded-atoms t
                                                     :assign-types t))
         (atom-table (chem:atom-table energy-function))
         (nonbond-db (chem:nonbond-force-field-for-aggregate atom-table))
         (ion1-topology (cando:lookup-topology ion1))
         (ion1-type (ion-topology-atom-type ion1-topology))
         (ion1-type-index (chem:find-type-index nonbond-db ion1-type))
         (ion1-ffnonbond (chem:get-ffnonbond-using-type-index nonbond-db ion1-type-index))
         (ion1-size (chem:get-radius-angstroms ion1-ffnonbond))
         (ion1-residue (chem:build-residue-single-name ion1-topology))
         (ion1-atom (chem:content-at ion1-residue 0))
         (ion1-mol (chem:make-molecule))
         (ion1-agg (chem:make-aggregate))
         (ion-vector (make-array (if ion2 (+ ion1-number ion2-number) ion1-number)))
         (ion-count 0)
         (ion-fail-counter 0)
         (target-charge 0.0)
         (ion2-size 0.0) 
         (ion2-type-index 0)
         ion2-ffnonbond ion2-topology ion2-residue ion2-atom ion2-mol ion2-agg)
    (chem:add-matter ion1-mol ion1-residue)
    (chem:add-matter ion1-agg ion1-mol)
    ;;Check the unit's validity
    (chem:map-atoms
     nil
     (lambda (r)
       (setf target-charge (+ (chem:get-charge r) target-charge)))
     aggregate)
    (if (= 0 target-charge)
        (progn
          (format t " ~s has a charge of 0. ~%" (chem:get-name aggregate))
          (if (= 0 ion1-number)
              (progn
                (error "Can't neutralize")
                (return-from add-ions-rand))
              (format t "Adding the ions anyway~%")))
        (format t "Charge ~,2f ~%" target-charge))
    ;;Make sure the ions are actually ions
    (if (= 0 (chem:get-charge ion1-atom))
        (progn
          (format t "~a is not an ion and is not appropriate for placement. ~%" ion1)
          (return-from add-ions-rand)))
    (if ion2
        (let* ((ion2-topology (cando:lookup-topology ion2))
               (ion2-type (ion-topology-atom-type ion2-topology)))
          (setf ion2-type-index (chem:find-type-index nonbond-db ion2-type))
          (setf ion2-ffnonbond (chem:get-ffnonbond-using-type-index nonbond-db ion2-type-index))
          (setf ion2-size (chem:get-radius-angstroms ion2-ffnonbond))
          (setf ion2-topology (cando:lookup-topology ion2))
          (setf ion2-residue (chem:build-residue-single-name ion2-topology))
          (setf ion2-atom (chem:content-at ion2-residue 0))
          (setf ion2-mol (chem:make-molecule))
          (setf ion2-agg (chem:make-aggregate))
          (chem:add-matter ion2-mol ion2-residue)
          (chem:add-matter ion2-agg ion2-mol)
          (if (= 0 (chem:get-charge ion2-atom))
              (error "~a is not an ion and is not appropriate for placement. ~%" ion2))))
    ;;Check validity of neutralization
    (if (= ion1-number 0)
        (progn
          (when (or (and (< (chem:get-charge ion1-atom) 0)
                       (< target-charge 0))
                  (and (> (chem:get-charge ion1-atom) 0)
                       (> target-charge 0)))
            (warn "1st Ion and target aggregate have charges of the same sign - aggregate charge: ~f; ion charge: ~f"
                  target-charge (chem:get-charge ion1-atom))
            (return-from add-ions-rand nil))
          ;;Get the nearest integer number of ions that we need to add to get as close as possible to neutral.
          (setf ion1-number (round (/ (abs target-charge) (abs (chem:get-charge ion1-atom)))))
          (when ion2
            (warn "Neutralization - can't do 2nd ion.")
            (return-from add-ions-rand nil))
          (format t "~d ~a ions required to neutraize. ~%" ion1-number (chem:get-name ion1-atom))))
    (if (not (chem:first-solvent-molecule-nspsol-bound-p atom-table))
        (error "No solvent present. Add solvent first.~%"))
    (let* ((final-solute-residue (chem:final-solute-residue-iptres atom-table))
           (total-number-of-molecules-nspm (chem:total-number-of-molecules-nspm atom-table))
           (first-solvent-molecule-nspsol (chem:first-solvent-molecule-nspsol atom-table))
           (total-number-of-solvent-molecules (+ (- total-number-of-molecules-nspm first-solvent-molecule-nspsol) 1)))
      (if ion2
          (if (>= (+ ion1-number ion2-number) total-number-of-solvent-molecules)
              (error "Too few solvent molecules to add ions.~%"))
          (if (>= ion1-number total-number-of-solvent-molecules)
              (error "Too few solvent molecules to add ions.~%")))
      (format t "Adding ~a counter ions to ~a. ~a solvent molecules will remain.~%"
              (if ion2 (+ ion1-number ion2-number) ion1-number) (chem:get-name aggregate)
              (if ion2 (- total-number-of-solvent-molecules ion1-number ion2-number)
                  (- total-number-of-solvent-molecules ion1-number)))
      ;;now actually add the ions
      (loop
        for ion1-copy = (chem:matter-copy ion1-agg)
        while (if ion2
                  (or (> ion1-number 0)
                      (> ion2-number 0))
                  (> ion1-number 0))
        do (if (> ion1-number 0)
               ;;Pick random solvent molecule to replace
               ;;Get position of solvent residue atom
               (let* ((random (random total-number-of-solvent-molecules))
                      (residue-count 0)
                      (picked-molecule-count (+ random first-solvent-molecule-nspsol))
                      (fail-counter 0)
                      (place-ion t)
                      picked-molecule picked-residue position)
                 (format t "picked-molecule-count ~a~%" picked-molecule-count)
                 (chem:map-molecules
                  nil
                  (lambda (m)
                    (when (= residue-count picked-molecule-count)
                      (setf picked-molecule m)
                      (setf picked-residue (chem:content-at m 0))
                      (setf position (chem:get-position (chem:content-at picked-residue 0))))
                    (incf residue-count))
                  aggregate)
                 (loop for i below ion-count
                       for point = (chem:get-position (chem:content-at (chem:content-at (chem:content-at (aref ion-vector i) 0) 0) 0))
                       for ion-distance = (geom:calculate-distance point position)
                       do (if (< ion-distance separation)
                              (progn
                                (setf place-ion nil)
                                (incf fail-counter))))
                 (if place-ion
                     (progn
                       (format t "~a: Placed ~a in ~a at (~,2f, ~,2f, ~,2f).~%"
                               ion-count ion1 (chem:get-name aggregate)
                               (geom:vx position) (geom:vy position) (geom:vz position))
                       ;;Save this ion's position if desired
                       ;;Copy ion unit, position, and add it to the unit
                       (chem:apply-transform-to-atoms ion1-copy (geom:make-m4-translate position))
                       (chem:add-matter aggregate ion1-copy)
                       (setf (aref ion-vector ion-count) ion1-copy)
                       (incf ion-count)
                       ;; Delete the solvent residue that was replaced
                       (chem:remove-atom-delete-bonds picked-residue (chem:content-at picked-residue 2))
                       (chem:remove-atom-dont-delete-bonds picked-residue (chem:content-at picked-residue 1))
                       (chem:remove-residue picked-molecule picked-residue)
                       (chem:remove-molecule aggregate picked-molecule)
                       (decf ion1-number)))
                 (if (> fail-counter 100)
                     (error "Impossible to place ~a ions with minumum separation of ~a~%"
                            (+ ion1-number ion2-number) separation))))
        do (if (and ion2 (> ion2-number 0))
               ;;Pick random solvent molecule to replace
               ;;Get position of solvent residue atom
               (let* ((ion2-copy (chem:matter-copy ion2-agg))
                      (random (random total-number-of-solvent-molecules))
                      (residue-count 0)
                      (picked-molecule-count (+ random first-solvent-molecule-nspsol))
                      (fail-count 0)
                      (place-ion t)
                      picked-molecule picked-residue position)
                 (format t "picked-molecule-count ~a~%" picked-molecule-count)
                 (chem:map-molecules
                  nil
                  (lambda (m)
                    (when (= residue-count picked-molecule-count)
                      (setf picked-molecule m)
                      (setf picked-residue (chem:content-at m 0))
                      (setf position (chem:get-position (chem:content-at picked-residue 0))))
                    (incf residue-count))
                  aggregate)          
                 (loop for i below ion-count
                       for point = (chem:get-position (chem:content-at (chem:content-at (chem:content-at (aref ion-vector i) 0) 0) 0))
                       for ion-distance = (geom:calculate-distance point position)
                       do (if (< ion-distance separation)
                              (progn
                                (setf place-ion nil)
                                (incf fail-count))))
                 (if place-ion
                     (progn
                       (format t "~a: Placed ~a in ~a at (~,2f, ~,2f, ~,2f).~%"
                               ion-count ion2 (chem:get-name aggregate)
                               (geom:vx position) (geom:vy position) (geom:vz position))
                       ;;Save this ion's position if desired
                       ;;Copy ion unit, position, and add it to the unit
                       (chem:apply-transform-to-atoms ion2-copy (geom:make-m4-translate position))
                       (chem:add-matter aggregate ion2-copy)
                       (setf (aref ion-vector ion-count) ion2-copy)
                       (incf ion-count)
                       ;; Delete the solvent residue that was replaced
                       (chem:remove-atom-delete-bonds picked-residue (chem:content-at picked-residue 2))
                       (chem:remove-atom-dont-delete-bonds picked-residue (chem:content-at picked-residue 1))
                       (chem:remove-residue picked-molecule picked-residue)
                       (chem:remove-molecule aggregate picked-molecule)
                       (decf ion2-number)))
                 (if (> fail-count 100)
                     (error "Impossible to place ~a ions with minumum separation of ~a~%"
                            (+ ion1-number ion2-number) separation)))))
             
        )))
