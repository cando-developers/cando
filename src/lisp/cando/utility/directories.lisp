(in-package :cando-utility)

(defun home/ (path)
  (if (member :docker *features*)
      (pathname (format nil "/work/~a" (namestring path)))
      (let ((env-home (ext:getenv "HOME")))
        (unless env-home (setf env-home "~"))
        (pathname (format nil "~a/~a" env-home (namestring path))))))


(defun current-directory ()
  "Return the current directory"
  (ext:getcwd))

(defun mkdir (dir)
  "Create a directory realtive to the users home directory"
  (ensure-directories-exist (merge-pathnames (uiop:ensure-directory-pathname dir))))

(defun calculate-directory (dir)
  (let ((pn-dir (pathname (uiop:ensure-directory-pathname dir))))
    (unless (eq :relative (car (pathname-directory pn-dir)))
      (error "Only relative directories are allowed"))
    (merge-pathnames pn-dir cando::*root-directory*)))
  
(defun set-current-directory (dir)
  "Set the current directory"
  (ext:chdir dir t))

(defun directory-files (&optional dir)
  (unless dir
    (setf dir (current-directory)))
  (let ((abs-dir (calculate-directory dir)))
    (uiop:directory-files abs-dir)))


