(in-package :cando-utility)

(defun mkdir (dir)
  (ensure-directories-exist (uiop:ensure-directory-pathname dir)))
  
