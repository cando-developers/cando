(in-package :cando-utility)

(defun current-directory ()
  "Return the current directory relative to the users home directory"
  (enough-namestring *default-pathname-defaults* cando::*root-directory*))

(defun mkdir (dir)
  "Create a directory realtive to the users home directory"
  (ensure-directories-exist (merge-pathnames (uiop:ensure-directory-pathname dir))))

(defun calculate-directory (dir)
  (let ((pn-dir (pathname (uiop:ensure-directory-pathname dir))))
    (unless (eq :relative (car (pathname-directory pn-dir)))
      (error "Only relative directories are allowed"))
    (merge-pathnames pn-dir cando::*root-directory*)))
  
(defun set-current-directory (dir)
  "Set the current directory (relative to the users home directory) to dir."
  (let ((abs-dir (calculate-directory dir)))
    (setf *default-pathname-defaults* abs-dir))
  (current-directory))

(defun directory-files (&optional dir)
  (unless dir
    (setf dir (current-directory)))
  (let ((abs-dir (calculate-directory dir)))
    (uiop:directory-files abs-dir)))


