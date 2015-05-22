


(defun mfix ()
  (interactive)
  (let ((cmd (concat "cat " (buffer-file-name) ".need")))
    (shell-command cmd 1)))