(in-package #:topology)

(defclass binned-shape-key-map-entry (cando.serialize:serializable)
  ((key :initarg :key :reader binned-map-entry-key)
   (value :initarg :value :accessor binned-map-entry-value)))

(defmethod initialize-instance :after ((entry binned-shape-key-map-entry) &key)
  (check-type (binned-map-entry-key entry) binned-shape-key)
  (setf (slot-value entry 'key)
        (make-binned-shape-key (binned-shape-key-entries (binned-map-entry-key entry)))))

(defclass binned-shape-key-map (cando.serialize:serializable)
  ((records :initarg :records :initform nil)
   ;; No initarg: this derived index is rebuilt, not serialized.
   (index :initform (make-hash-table :test #'equal)))
  (:documentation "Content-indexed shape keys. Stored keys are private snapshots.
The hash table's structural keys are an internal detail, not the public API."))

(defmethod initialize-instance :after ((map binned-shape-key-map) &key)
  (unless (binned-shape-key-proper-list-p (slot-value map 'records))
    (error "Shape-key map records must be a proper list"))
  (let ((index (slot-value map 'index)))
    (dolist (record (slot-value map 'records))
      (check-type record binned-shape-key-map-entry)
      (let ((entries (binned-shape-key-entries (binned-map-entry-key record))))
        (when (gethash entries index)
          (error "Duplicate key in serialized shape-key map: ~s" entries))
        (setf (gethash entries index) record)))))

(defun make-binned-shape-key-map ()
  (make-instance 'binned-shape-key-map))

(defun binned-shape-key-map-keys (map)
  "Return independent key objects for inspection."
  (let (keys)
    (map-binned-shape-keys (lambda (key value) (declare (ignore value)) (push key keys)) map)
    keys))

(defun binned-shape-key-map-value (key map &optional default)
  "Return VALUE and PRESENT-P for an equivalent key, as with GETHASH."
  (check-type key binned-shape-key)
  (multiple-value-bind (record presentp)
      (gethash (binned-shape-key-entries key) (slot-value map 'index))
    (values (if presentp (binned-map-entry-value record) default) presentp)))

(defun (setf binned-shape-key-map-value) (value key map &optional default)
  (declare (ignore default))
  (check-type key binned-shape-key)
  (let* ((entries (binned-shape-key-entries key))
         (index (slot-value map 'index))
         (record (gethash entries index)))
    (if record
        (setf (binned-map-entry-value record) value)
        (let ((new (make-instance 'binned-shape-key-map-entry :key key :value value)))
          (push new (slot-value map 'records))
          (setf (gethash entries index) new))))
  value)

(defun binned-shape-key-map-count (map)
  (hash-table-count (slot-value map 'index)))

(defun map-binned-shape-keys (function map)
  "Call FUNCTION with each selected-key object and its value.
Keys are copies, so modifying a callback key cannot corrupt the map."
  (dolist (record (slot-value map 'records))
    (funcall function
             (make-binned-shape-key (binned-shape-key-entries (binned-map-entry-key record)))
             (binned-map-entry-value record)))
  nil)
