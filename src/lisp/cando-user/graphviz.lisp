(in-package :cando-user)

(defun graphviz-dot (graphviz-renderer-name object-to-render &optional filename)
  "Generate a PDF for an object using graphviz dot.  
The output will be written to a file named **filename**.  
If **filename** is nil then it will be written to the name graphviz-renderer-name.pdf
Graphviz dot must be installed and in the current path."
  (let* ((pdf-pname (if filename
                        (make-pathname :name (pathname-name filename) :type "pdf")
                        (pathname (format nil "~a.pdf" (string graphviz-renderer-name)))))
         (dot-pname (make-pathname :type "dot" :defaults pdf-pname))
         (dot-cmd (list "dot" "-o" (namestring pdf-pname)
                            "-T" "pdf"
                            (namestring dot-pname))))
    (format t "Running:  ~{~a ~}~%" dot-cmd)
    (funcall graphviz-renderer-name object-to-render dot-pname)
    (ext:vfork-execvp dot-cmd)))
