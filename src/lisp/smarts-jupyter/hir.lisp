(in-package :smarts-jupyter)

(defun hir-graph (hir)
  (format t "Drawing graph ~a~%" hir)
  (let ((nodes nil)
        (root-id nil)
        (instruction-to-node (make-hash-table))
        (instruction-to-node-id (make-hash-table)))
    (cleavir-ir:map-instructions
     (lambda (instruction)
       (let* ((node-id (gensym))
              (label (class-name (class-of instruction)))
              (node (make-instance 'cytoscape:element
                                   :group "nodes"
                                   :data (j:make-object "id" (string node-id)
                                                        "label" label))))
         (push node nodes)
         (unless root-id
           (setf root-id (string node-id)))
         (setf (gethash instruction instruction-to-node) node
               (gethash instruction instruction-to-node-id) node-id)))
     hir)
    (let ((edges nil))
      (cleavir-ir:map-instructions
       (lambda (instruction)
         (loop for succ in (cleavir-ir:successors instruction)
               for edge = (make-instance 'cytoscape:element
                                         :group "edges"
                                         :data (j:make-object "source" (string (gethash instruction instruction-to-node-id))
                                                              "target" (string (gethash succ instruction-to-node-id))))
               do (push edge edges)))
       hir)
      (format t "Number of nodes: ~a  edges: ~a~%" (length nodes) (length edges))
      (make-instance 'cytoscape:cytoscape-widget
                     :graph-layouts (list (make-instance 'cytoscape:dagre-layout
                                                         :roots (list root-id)))
                     :graph-style "* { label: data(label); font-size: 10; }"
                     :elements (append (nreverse nodes) (nreverse edges))
                     :layout (make-instance 'jupyter-widgets:layout :width "auto" :height "2000px")))))


(defun draw-form-hir2 (form)
  "Generate a HIR graph for the form using the cst compiler"
  (let (result)
    (cmp::with-compiler-env ()
      (let* ((module (cmp::create-run-time-module-for-compile)))
        ;; Link the C++ intrinsics into the module
        (let* ((cleavir-cst-to-ast:*compiler* 'cl:compile)
               (cst (cst:cst-from-expression form))
               (ast (cleavir-cst-to-ast:cst-to-ast cst nil clasp-cleavir::*clasp-system*))
               (hoisted-ast (clasp-cleavir::hoist-ast ast))
               (hir (clasp-cleavir::ast->hir hoisted-ast)))
          (setq result hir)
          (cmp::llvm-create-module "foo"))))
    (hir-graph result)))


(defun draw-form-hir (form)
  "Generate a HIR graph for the form using the cst compiler"
  (let ((clasp-cleavir::*save-hir* t))
    (compile 'nil form)
    (hir-graph clasp-cleavir::*hir*)))
