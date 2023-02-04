(in-package :cando-user)

(defun sharp-!-reader (stream subchar arg)
  "Read to the end of a line honoring continuation characters \\ and accumulate a string
   that is passed to a bash subshell.  The output from the evaluation of the command
   in the subshell is printed to the standard-output.
   Note: When a change directory command like #!cd <xxx> is evaluated, when the subshell returns
   the current working directory will not be changed."
  (declare (ignore arg subchar))
  (let ((accumulated (loop named accumulate
                           with sout = (make-string-output-stream)
                           for line = (read-line stream)
                           for trim-line = (string-right-trim (list #\space #\tab) line)
                           for len = (length trim-line)
                           for last-char = (if (= len 0) nil (elt trim-line (1- len)))
                           if (and (char= last-char #\\))
                             do (princ (subseq trim-line 0 (1- len)) sout)
                           else
                             do (princ trim-line sout)
                           unless (and last-char (char= last-char #\\))
                             do (return-from accumulate (get-output-stream-string sout))
                           do (princ #\space sout))))
    (multiple-value-bind (instream return-code process)
        (ext:run-program "/bin/bash" (list "-c" accumulated))
      (declare (ignore return-code process))
      (loop named print-output
            for line = (read-line instream nil nil)
            if line
              do (format t "~a~%" line)
            else
              do (return-from print-output nil)))
    (values)))

(eval-when (:execute :load-toplevel)
  ;; Install the dispatching reader macro
  (set-dispatch-macro-character #\# #\! 'sharp-!-reader))

