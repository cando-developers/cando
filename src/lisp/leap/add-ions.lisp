(in-package :leap.add-ions)

(defun add-ions (mol ion1-element ion1-number &optional ion2-element ion2-number)
  (let* ((ion1 (chem:make-atom ion1-element ion1-element))
         (res-ion1 (chem:make-residue ion1-element))
         (mol-ion1 (chem:make-molecule))
         (agg-ion1 (chem:make-aggregate))
         (target-charge 0)
         (atom-radius 1.5)
         (ion1-size 1.37) ;;temporally 
         (ion2-size 0.0) 
         (ion-min-size 0.0)
         (oct-shell 1)
         (grid-space 1.0)
         (shell-extent 4.0)
         (at-octree 1)
         (dielectric 1))
    (chem:set-charge ion1 '1.0) ;;temporally
    (chem:add-matter res-ion1 ion1)
    (chem:add-matter mol-ion1 res-ion1)
    (chem:add-matter agg-ion1 mol-ion1)
    (if (and ion2-element (= 0 ion2-number))
        (error "'0' is not allowed as the value for the second ion."))
    ;;Consider target unit's charge
    (chem:map-atoms
     nil
     (lambda (r)
       (setf target-charge (+ (chem:get-charge r) target-charge)))
     mol)
    (if (and (= 0 target-charge)
             (= 0 ion1-number))
        (error "Can't neutralize.")
        (format t "Total charge ~f~%" target-charge))
    ;;Consider neutralizetion    
    (if (= ion1-number 0)
        (progn
          (if (or (and (< (chem:get-charge ion1) 0)
                       (< target-charge 0))
                  (and (> (chem:get-charge ion1) 0)
                     (> target-charge 0)))
              (error "1st ion & target are same charge"))
          ;;Get the nearest integer number of ions that we need to add to get as close as possible to neutral.
          (setf ion1-number (round (/ (abs target-charge) (chem:get-charge ion1))))
          (if ion2-element
              (error "Neutralization - can't do 2nd ion."))))
    (format t "~d ~a ions required to neutraize. ~%" ion1-number (chem:get-name ion1))
    ;;Consider ion sizes and postions
    (if ion2-element
        (setf ion-min-size (min ion1-size ion2-size))
        (setf ion-min-size ion1-size))
    ;;Build grid and calc potential on it
    (let ((oct-tree (core:make-cxx-object 'chem:octree)))
      (chem:oct-oct-tree-create oct-tree mol oct-shell grid-space ion-min-size shell-extent 0 t)
      (multiple-value-bind (min-charge-point max-charge-point)
          (chem:oct-tree-init-charges oct-tree at-octree dielectric ion1-size)
        (loop for i from 0 below ion1-number
           for new-point = (if (< (chem:get-charge ion1) 0)
                               (setf new-point max-charge-point)
                               (setf new-point min-charge-point))
           for ion1-copy = (chem:matter-copy agg-ion1)
           for ion1-transform = (geom:make-m4-translate new-point)
           do (chem:apply-transform-to-atoms ion1-copy ion1-transform)
           do (chem:add-matter mol ion1-copy)
           do (format t "Placed ~a in ~a at ~a ~a ~a~%" (chem:get-name ion1) 
                      (chem:get-name mol)
                      (geom:vx new-point)
                      (geom:vy new-point)
                      (geom:vz new-point))
           do (chem:oct-tree-delete-sphere oct-tree new-point (if ion2-element
                                                                                  (+ ion1-size ion2-size)
                                                                                  (+ ion1-size ion1-size)))
           do (multiple-value-bind (min-new-charge-point max-new-charge-point)
                  (chem:oct-tree-update-charge oct-tree new-point (chem:get-charge ion1)
                                               (if ion2-element ion2-size ion1-size))
                (setf min-charge-point min-new-charge-point)
                (setf max-charge-point max-new-charge-point))
             )
       #+(or) (if ion2-element
            (loop for i from 0 below ion2-number
               for new-point = (if (< (chem:get-charge ion2) 0)
                                   (setf new-point max-charge-point)
                                   (setf new-point min-charge-point))
               for ion2-copy = (chem:matter-copy agg-ion2)
               for ion2-transform = (geom:make-m4-translate new-point)
               do (chem:apply-transform-to-atoms ion2-copy ion2-transform)
               do (chem:add-matter mol ion2-copy)
               do (format t "Placed ~a in ~a at ~a ~a ~a~%" (chem:get-name ion2) 
                          (chem:get-name mol)
                          (geom:vx new-point)
                          (geom:vy new-point)
                          (geom:vz new-point))
               do (chem:oct-tree-delete-sphere oct-tree oct-tree-solute new-point (if ion1-element
                                                                                      (+ ion2-size ion1-size)
                                                                                      (+ ion2-size ion2-size)))
               do (multiple-value-bind (min-new-charge-point max-new-charge-point)
                      (chem:oct-tree-update-charge oct-tree oct-tree-solute new-point (chem:get-charge ion2)
                                                   (if ion1-element ion1-size ion2-size))
                    (setf min-charge-point min-new-charge-point)
                    (setf max-charge-point max-new-charge-point))))
        ))
    
    )
  )

