(in-package :design.topology)


(defun new-joint-template-factory (parent-template node in-plug)
  (let* ((atom-name (topology:name node))
         (child-indexes (mapcar #'topology:constitution-atom-index (topology:children node)))
         (constitution-atoms-index (topology:constitution-atom-index node))
         (gparent-template (if parent-template
                               (design.joint-tree:parent parent-template)
                               nil))
         (ggparent-template (if gparent-template
                                (design.joint-tree:parent gparent-template)
                                nil)))
    (cond
      ((and (null parent-template) (typep in-plug 'chem:in-plug))
       (design.joint-tree:make-in-plug-bonded-joint-template constitution-atoms-index
                                           :atom-name atom-name
                                           :parent nil
                                           :in-plug in-plug))
      ((typep in-plug 'chem:in-plug)
       (design.joint-tree:make-bonded-joint-template constitution-atoms-index
                                   :atom-name atom-name
                                   :parent parent-template))
      ((null parent-template)
       (design.joint-tree:make-jump-joint-template constitution-atoms-index
                                 :atom-name atom-name
                                 ))
      ((null gparent-template)
       (let ((stub-joints (cond
                            ((eql 0 (first child-indexes))
                             (list parent-template))
                            ((eql 1 (first child-indexes))
                             (list parent-template (design.joint-tree:sibling parent-template 0)))
                            (t
                             (list parent-template
                                   (design.joint-tree:sibling parent-template 1)
                                   (design.joint-tree:sibling parent-template 0))))))
         (design.joint-tree:make-complex-bonded-joint-template constitution-atoms-index
                                             :atom-name atom-name
                                             :stub-joints stub-joints)))
      ((null ggparent-template)
       (let ((stub-joints (cond
                            ((and (eql 0 (first child-indexes)) (eql 0 (second child-indexes)))
                             (list parent-template
                                   gparent-template))
                            ((and (> 0 (first child-indexes)) (eql 0 (second child-indexes)))
                             (list parent-template
                                   gparent-template
                                   (design.joint-tree:sibling parent-template 0)))
                            (t
                             (list parent-template
                                   gparent-template
                                   (design.joint-tree:sibling gparent-template 0))))))
         (design.joint-tree:make-complex-bonded-joint-template constitution-atoms-index
                                                               :atom-name atom-name
                                                               :stub-joints stub-joints)))
      (t (design.joint-tree:make-bonded-joint-template constitution-atoms-index
                                                       :atom-name atom-name
                                                       :parent parent-template)))))


(defun build-joint-template-recursively (parent root in-plug constitution-atoms constitution-name topology-name )
  (let ((root-template (new-joint-template-factory parent root in-plug constitution-atoms constitution-name topology-name))
        (children (topology:children root)))
    (format t "root-template = ~s~%" root-template)
    (loop for child in children
          for sub-child-index from 0
          do (let ((child-template (build-atom-tree-template-recursively root-template
                                                                         child
                                                                         in-plug
                                                                         constitution-atoms
                                                                         constitution-name
                                                                         topology-name )))
               (design.joint-tree:add-child root-template child-template)))
    root-template))



(defun build-joint-template (graph)
  (let ((root-node (topology:root-node graph)))
    (format t "root-node = ~a~%" root-node)
    (build-joint-template-recursively nil root-node)
    nil))



