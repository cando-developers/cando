;;;; lsp-builder-generator.el --- TODO


;;; Commentary:

(require 'derived)
(require 'lsp-mode)
(require 'lsp-ui)


;;; Code:

;;; Language Server Protocol

(lsp-define-stdio-client
 smarts nil
 #'(lambda () default-directory)
 '("/home/jmoringe/code/cl/language.smarts/language-server")
 :language-id-fn (lambda (buffer)
                   (with-current-buffer buffer
                     (ecase major-mode
                       (smiles-mode "smiles")
                       (smarts-mode "smarts"))))
 :docstring
 "Language Server for SMILES and SMARTS languages.")

;;; Language Modes

(defmacro define-smarts-language-mode
    (function-name name file-extension)
  (let ((hook-name      (intern (concat (symbol-name function-name) "-hook")))
        (filename-regex (concat "\\." file-extension "\\'")))
    `(progn
       (define-derived-mode ,function-name yaml-mode ,name)
       (add-to-list 'auto-mode-alist '(,filename-regex . ,function-name))
       (add-hook ',hook-name
                 (lambda ()
                   (require 'company-capf)
                   (setf company-backends              '((company-capf))
                         company-minimum-prefix-length 1)
                   (smarts-enable)))
       (lsp-ui-flycheck-add-mode ',function-name))))

(define-smarts-language-mode smarts-mode
  "Smarts" "smarts")

(define-smarts-language-mode smiles-mode
  "Smiles" "smiles")

;;; test.el ends here
