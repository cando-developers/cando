(in-package :cando-utility)


(defun mkdir (dir)
  (ensure-directories-exist (merge-pathnames (uiop:ensure-directory-pathname dir))))

(defun set-current-directory (dir)
  (setf *default-pathname-defaults* (uiop:ensure-directory-pathname dir)))
  
