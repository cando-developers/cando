;;; Lexical stuff

(macrolet ((define-delimiter (character &optional (name (string character)))
             (let ((rule-name (symbolicate '#:delimiter- name)))
               `(defrule ,rule-name
                    ,character
                  (:constant nil)))))
  (define-delimiter #\()
  (define-delimiter #\)))

;;; Actual grammar
