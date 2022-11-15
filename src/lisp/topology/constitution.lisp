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

(defclass topology ()
  ((name :initarg :name :accessor name)
   (constitution :initarg :constitution :accessor constitution)
   (property-list :initform nil :initarg :property-list :accessor property-list)
   (plugs :initarg :plugs :accessor plugs)
   (stereoisomer-atoms :initarg :stereoisomer-atoms :accessor stereoisomer-atoms)
   (default-out-plug-name :initarg :default-out-plug-name :accessor default-out-plug-name)))

(defmethod print-object ((obj topology) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~a" (name obj))))
