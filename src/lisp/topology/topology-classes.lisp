(in-package :topology)

(defclass constitution-atom ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (index :initarg :index :accessor index)
   (element :initarg :element :accessor element)
   (atom-type :initarg :atom-type :accessor atom-type)
   (stereochemistry-type :initarg :stereochemistry-type :accessor stereochemistry-type)
   (properties :initform nil :initarg :properties :accessor properties)
   (bonds :initform nil :initarg :bonds :accessor bonds)))

(defmethod print-object ((obj constitution-atom) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (atom-name obj))))

(defclass constitution-bond ()
  ((to-atom-index :initarg :to-atom-index :accessor to-atom-index)
   (order :initarg :order :accessor order)))

(defclass constitution-atoms ()
  ((atoms :initarg :atoms :accessor atoms)))

(defclass stereoconfiguration ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (configuration :initarg :configuration :accessor configuration)))

(defmethod print-object ((obj stereoconfiguration) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (atom-name obj) (configuration obj))))

(defclass stereoisomer ()
  ((name :initarg :name :accessor name)
   (pdb :initarg :pdb :accessor pdb)
   (stereoisomer-index :initarg :stereoisomer-index :accessor stereoisomer-index)
   (configurations :initarg :configurations :accessor configurations)))

(defmethod print-object ((obj stereoisomer) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~s" (name obj) (configurations obj))))

(defclass constitution ()
  ((name :initarg :name :accessor name)
   (constitution-atoms :initarg :constitution-atoms :accessor constitution-atoms)
   (plugs :initarg :plugs :accessor plugs)
   (topology-list :initform nil :initarg :topology-list :accessor topology-list)
   (stereo-information :initarg :stereo-information :accessor stereo-information)))

(defclass stereoisomer-atom ()
  ((atom-name :initarg :atom-name :accessor atom-name)
   (constitution-atom-index :initarg :constitution-atom-index :accessor constitution-atom-index)
   (atom-charge :initarg :atom-charge :accessor atom-charge)
   (atom-type :initarg :atom-type :accessor atom-type)))

(defclass stereoisomer-virtual-atom (stereoisomer-atom)
  ())

(defclass topology ()
  ((name :initarg :name :accessor name)
   (constitution :initarg :constitution :accessor constitution)
   (property-list :initform nil :initarg :property-list :accessor property-list)
   (plugs :initarg :plugs :accessor plugs)
   (joint-template :initarg :joint-template :accessor joint-template)
   (stereoisomer-atoms :initform (make-array 4 :adjustable t :fill-pointer 0)
                       :initarg :stereoisomer-atoms :accessor stereoisomer-atoms)
   ))

(defmethod print-object ((obj topology) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (name obj))))

(defun has-plug-named (topology plug-name)
  (gethash plug-name (plugs topology)))

(defun plug-named (topology plug-name)
  (gethash plug-name (plugs topology)))

(defun find-in-plug (topology)
  (maphash (lambda (name plug)
             (declare (ignore name))
             (when (typep plug 'in-plug)
               (return-from find-in-plug plug)))
           (plugs topology)))

(defun out-plugs-as-list (topology)
  (let (out-plugs)
    (maphash (lambda (name plug)
               (declare (ignore name))
               (when (typep plug 'out-plug)
                 (push plug out-plugs)))
             (plugs topology))
    out-plugs))

(defun all-out-plug-names-that-match-in-plug-name (topology in-plug-name)
  (let (out-plugs
        (in-coupling-name (coupling-name in-plug-name)))
    (maphash (lambda (name plug)
               (declare (ignore name))
               (when (and (typep plug 'out-plug)
                          (eq (coupling-name (name plug))
                              in-coupling-name))
                 (push (name plug) out-plugs)))
             (plugs topology))
    out-plugs))

(defun plugs-as-list (topology)
  (let (plugs)
    (maphash (lambda (name plug)
               (declare (ignore name))
               (push plug plugs))
             (plugs topology))
    plugs))

(defclass plug ()
  ((name :initarg :name :accessor name)
   (atom-names :initform (make-array 16)
               :accessor atom-names)))

(defclass in-plug (plug)
  ())

(defclass out-plug (plug)
  ())

(defmethod print-object ((obj plug) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a" (name obj) (atom-names obj))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Monomer, coupling, and oligomer
;;;


(defclass monomer ()
  ((id :initarg :id :accessor id)
   (couplings :initform (make-hash-table)
              :accessor couplings)
   (current-stereoisomer-offset :initform 0 :accessor current-stereoisomer-offset)
   (monomers :initarg :monomers :accessor monomers)))

(defmethod print-object ((obj monomer) stream)
   (print-unreadable-object (obj stream :type t)
    (format stream "~a" (monomers obj))))

(defun number-of-stereoisomers (monomer)
  (length (monomers monomer)))

(defun current-stereoisomer-name (monomer)
  (elt (monomers monomer) (current-stereoisomer-offset monomer)))

(defun get-current-topology (monomer)
  (let ((name (elt (monomers monomer) (current-stereoisomer-offset monomer))))
    (chem:find-topology name t)))

(defun has-in-coupling-p (monomer)
  (maphash (lambda (plug-name plug)
             (declare (ignore plug-name))
             (when (typep plug 'topology:in-plug)
               (return-from has-in-coupling-p t)))
           (couplings monomer)))

(defun monomer-plug-named (monomer plug-name)
  (gethash plug-name (couplings monomer)))

(defclass coupling ()
  ())

(defclass directional-coupling (coupling)
  (
   (source-plug-name :initarg :source-plug-name :accessor source-plug-name)
   (target-plug-name :initarg :target-plug-name :accessor target-plug-name)
   (source-monomer :initarg :source-monomer :accessor source-monomer)
   (target-monomer :initarg :target-monomer :accessor target-monomer)))

(defmethod print-object ((obj directional-coupling) stream)
   (print-unreadable-object (obj stream :type t)
    (format stream "~a ~a ~a ~a" (source-monomer obj) (source-plug-name obj) (target-plug-name obj) (target-monomer obj))))


(defclass ring-coupling (coupling)
  (
   (plug1 :initarg :plug1 :accessor plug1)
   (plug2 :initarg :plug2 :accessor plug2)
   (monomer1 :initarg :monomer1 :accessor monomer1)
   (monomer2 :initarg :monomer2 :accessor monomer2)))

(defclass oligomer ()
  ((monomers :initform (make-array 16 :adjustable t :fill-pointer 0)
             :initarg :monomers :accessor monomers)
   (couplings :initform (make-array 16 :adjustable t :fill-pointer 0)
              :initarg :couplings :accessor couplings)))

(defun add-monomer (oligomer monomer)
  (vector-push-extend monomer (monomers oligomer)))

(defun number-of-sequences (oligomer)
  (let ((num 1))
  (loop for monomer across (monomers oligomer)
        do (setf num (* num (length (monomers monomer)))))
    num))

(defun goto-sequence (oligomer index)
  (let* ((bases (loop for monomer across (monomers oligomer)
                      collect (length (monomers monomer))))
         (digits (sys:positive-integer-to-mixed-base-digits index bases)))
    (loop for monomer across (monomers oligomer)
          for digit in digits
          do (setf (current-stereoisomer-offset monomer) digit))))
