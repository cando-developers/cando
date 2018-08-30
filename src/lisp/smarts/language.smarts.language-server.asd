(in-package :asdf-user)

(defsystem language.smarts.language-server
  :description "TODO"
  :license     "TODO"
  :author      "Jan Moringen <jmoringe@techfak.uni-bielefeld.de>"

  :version     (:read-file-form "version-string.sexp")
  :depends-on  ((:version "protocol.language-server"                     "0.1")

                (:version "language.smiles"                              (:read-file-form "version-string.sexp"))
                (:version "language.smarts"                              (:read-file-form "version-string.sexp"))

                (:version "text.source-location.lookup"                  "0.1")
                (:version "text.source-location.source-tracking-builder" "0.1"))

  :components  ((:module     "language-server"
                 :pathname   "src/language-server"
                 :serial     t
                 :components ((:file       "package"))))

  :entry-point "LANGUAGE.SMARTS.LANGUAGE-SERVER::LANGUAGE-SERVER")
