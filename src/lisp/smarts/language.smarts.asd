(in-package :asdf-user)

(defsystem "language.smarts"
  :version    (:read-file-form "version-string.sexp")
  :depends-on ("alexandria"

               (:version "esrap"                         "0.16")

               (:version "parser.common-rules"           "0.4")
               (:version "parser.common-rules.operators" "0.4")

               (:version "architecture.builder-protocol" "0.9")

               (:version "language.smiles"               (:read-file-form "version-string.sexp")))

  :components ((:module    "parser"
                :pathname  "src/smarts/parser"
                :serial    t
                :components ((:file       "package")
                             (:file       "grammar")
                             (:file       "protocol"))))

  :in-order-to ((test-op (test-op "language.smarts/test"))))

(defsystem "language.smarts/playground"
  :depends-on ("language.smarts"

               "architecture.builder-protocol.print-tree"))

(defsystem "language.smarts/test"
  :depends-on ("alexandria"

               "fiveam"

               "language.smarts")

  :components ((:module     "test"
                :pathname   "test"
                :serial     t
                :components ((:file       "package")))

               (:module     "parser"
                :depends-on ("test")
                :pathname   "test/smarts/parser"
                :serial     t
                :components ((:file       "package")
                             (:file       "grammar" ))))

  :perform    (test-op (operation component)
                (uiop:symbol-call '#:language.smarts.test '#:run-tests)))
