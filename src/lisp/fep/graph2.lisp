;; <script type=\"text/javascript\" src=\"https://www.graphdracula.net/js/graph.js\"></script>

;;; Copyright (c) 2019, Christian E. Schafmeister
;;; Published under the GPL 2.0.  See COPYING
;;;

(in-package :graph)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;
;;; Drawing a node in the stream

(defun render-node (stream id name &key (color "white") (shape "rectangle"))
  (format stream "  \"data\" : {~%")
  (format stream "    \"id\" : \"~a\"~%" id)
  (format stream "    label ~s~%" (format nil "~a" name))
  (format stream "    graphics~%")
  (format stream "    [~%")
  (format stream "      type ~s~%" (string-downcase (string shape)))
  (format stream "      fill ~s~%" (string-downcase (string color)))
  (format stream "    ]~%")
  (format stream "  ]~%"))

(defun render-morph (stream source target label &key type)
  (when source
    (format stream "  morph~%")
    (format stream "  [~%")
    (format stream "    source ~a~%" source)
    (format stream "    target ~a~%" target)
    (format stream "    label ~s~%" label)
    (format stream "    graphics~%")
    (format stream "    [~%")
    (format stream "      width 2~%")
    (format stream "      type ~s~%" "line")
    (format stream "    ]~%")
    (format stream "  ]~%")))

(defun render-node-json (stream id name &key (color "white") (shape "rectangle") last-node)
  (format stream "  {~%")
  (format stream "  \"data\" : {~%")
  (format stream "    \"id\" : \"~a\",~%" id)
  (format stream "    \"label\" : ~s,~%" (format nil "~a" name))
  (format stream "    \"graphics\" : {~%")
  (format stream "      \"type\" : ~s,~%" (string-downcase (string shape)))
  (format stream "      \"fill\" : ~s~%" (string-downcase (string color)))
  (format stream "    }~%")
  (format stream "  }~%")
  (format stream "  }")
  (unless last-node (format stream ","))
  (terpri stream))

(defun render-morph-json (stream source target label &key (width 2) type last-morph)
  (when source
    (format stream "  {~%")
    (format stream "  \"data\" : {~%")
;;    (format stream "    \"id\" : \"~a\",~%" id)
    (format stream "    \"source\" : ~s,~%" (format nil "~a" source))
    (format stream "    \"target\" : ~s,~%" (format nil "~a" target))
    (format stream "    \"label\" : ~s,~%" (format nil "~a" label))
    (format stream "    \"graphics\" : {~%")
    (format stream "      \"width\" : ~s,~%" width)
    (format stream "      \"type\" : ~s~%" "line")
    (format stream "    }~%")
    (format stream "  }~%")
    (format stream "  }")
    (unless last-morph (format stream ","))
    (terpri stream)))



(defun generate-graph (jobs)
  (let ((nodes (make-hash-table)))
    (loop for node in (fep:nodes jobs)
          do (setf (gethash node nodes) (fep:name node)))
    (let ((node-string (with-output-to-string (sout)
                         (maphash (lambda (node name)
                                    (format sout "{data: { id: '~a', label: '~a' }},~%" (string name) (string name)))
                                  nodes)))
          (morph-string (with-output-to-string (sout)
                         (loop for (from to) in (fep:morphs jobs)
                               do (format sout "{data: { source: '~a', target: '~a' }},~%" (string (fep:name from)) (string (fep:name to)))))))
      (values node-string morph-string))))

    
(defparameter *graph-form*
"<div class=\"custom-widget\">
  <div class=\"canvas\" style=\"width: 100%; min-height: 500px; border: 1px solid gray;\"></div>
  <script type=\"text/javascript\">
    var widgets = document.querySelectorAll(\".custom-widget\");
    ((self)=>{
        require.config({
            paths: {
                cytoscape: \"https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.2.8/cytoscape.min\"
            }
        });
        require([\"cytoscape\"], (cytoscape)=>{
            var canvas = self.querySelector(\".canvas\");
            var cy = cytoscape({
                container: canvas,
                elements: {
                    nodes: [ ~a
                    ],
                    morphs: [ ~a
                    ]},
                
                style: [
                    {
                        selector: 'node',
                        style: {
                            'background-color': 'white',
                            'border-width': '1px',
                            'border-style': 'solid',
                            'border-color': '#ccc',
                            'padding': '5px',
                            'label': 'data(label)',
                            'text-valign': \"center\",
                            'width': \"label\"
                        }
                    },
                    {
                        selector: 'node:selected',
                        style: {
                            'background-color': '#ccc',
                        }
                    },
                    {
                        selector: 'morph',
                        style: {
                            'label': 'data(label)',
                            'curve-style': 'bezier',
                            'width': 1,
                            'line-color': '#ccc',
                            'target-arrow-shape': 'triangle',
                            'target-arrow-color': '#ccc'
                        }
                    }
                ],

                layout: {
                    name: 'cose',
                    animate: false
                }
            });
            // Add double tap event
            var tappedBefore;
            var tappedTimeout;
            cy.on('tap', (ev)=>{
                var tappedNow = ev.target;                
                if (tappedTimeout && tappedBefore) {
                    clearTimeout(tappedTimeout);
                }
                if(tappedBefore === tappedNow) {
                    tappedNow.trigger('doubleTap', [ev.position]);
                    tappedBefore = null;
                } else {
                    tappedTimeout = setTimeout(()=>{ tappedBefore = null; }, 300);
                    tappedBefore = tappedNow;
                }
            });
            // Event Handling
            canvas.querySelector(\"canvas:last-child\").setAttribute(\"tabindex\", 1);
            var shift = false;
            canvas.querySelector(\"canvas:last-child\").addEventListener(\"keydown\", (ev)=>{
                if(ev.shiftKey || ev.keyCode == 16) shift = true;
            });
            canvas.querySelector(\"canvas:last-child\").addEventListener(\"keyup\", (ev)=>{
                if(ev.shiftKey || ev.keyCode == 16) shift = false;
                switch(ev.keyCode){
                    case 46:
                        var selected = cy.$(':selected');
                        if(selected) selected.remove();
                        break;
                }
            });
            cy.on(\"tap\", (ev)=>{
                canvas.querySelector(\"canvas:last-child\").focus();
            });
            var shiftSelected = null;
            cy.on(\"tap\", \"node\", (ev)=>{
                if(shift){
                    if(shiftSelected != null){
                        console.log(shiftSelected, ev.target);
                        cy.add({group: \"edges\", data: {source: shiftSelected.id(), target: ev.target.id()}});
                    }
                }
                shiftSelected = ev.target;
            });
            cy.on(\"doubleTap\", (ev, position)=>{
                if(ev.target == cy){
                    var id = prompt(\"Enter the node ID\");
                    if(id) cy.add({group: \"nodes\", data: {id: id, label: id}, position: position});
                }else{
                    var label = prompt(\"Enter the new label\", ev.target.data(\"label\"));
                    if(label) ev.target.data(\"label\", label);
                }
            });
            canvas.querySelector(\"canvas:last-child\").focus();
            self.querySelector(\".save\").addEventListener(\"click\", () => {
             var command = '(defparameter graph::*graph* \"' + JSON.stringify(cy.json()).replace(/\"/g,'\\\\\"') + '\")';
             console.log(\"Executing Command: \" + command);
             var kernel = IPython.notebook.kernel;
             kernel.execute(command);
            });
        });
    })(widgets[widgets.length-1]);
  </script>
  <p>
    Double click an empty space to add a node.
    Double click a node or morph to change its label.
    Click a node, then shift-click another to add an morph.
    Select a node and hit Delete to remove it.
  </p>
<button class=\"save\">Save</button>
</div>")


#+(or)
(defun fep-graph (pairs)
  (multiple-value-bind (node-string morph-string)
    (generate-graph pairs)
    (let ((html-string
            (with-output-to-string (sout)
              (format sout *graph-form* node-string morph-string))))
      ;;;(format t "~a~%" html-string)
      (cl-jupyter-user:html html-string))))



(defun save-graph (graph &optional (filename "graph.cyjs"))
  (with-open-file (fout filename :direction :output)
    (format fout "{~%")
    (format fout "  \"elements\" : {~%")
    (format fout "    \"nodes\" : [~%")
    (loop for index from 0
          for last-node = (= index (1- (length (fep:nodes (fep:jobs graph)))))
          for node in (fep:nodes (fep:jobs graph))
          do (render-node-json fout (fep:name node) (fep:name node) :shape "circle" :last-node last-node))
    (format fout "  ],~%")
    (format fout "  \"edges\" : [~%")
    (loop for index from 0
          for last-morph = (= index (1- (length (fep:morphs (fep:jobs graph)))))
          for morph in (fep:morphs (fep:jobs graph))
          do (render-morph-json fout (fep:name (fep:source morph)) (fep:name (fep:target morph)) "dummy" :last-morph last-morph))
    (format fout "]~%")
    (format fout "}~%")
    (format fout "}~%")
    )
  graph)

#||
;Test with

(defparameter *pairs* '((:aaba :aabb) (:aabb :aabc) (:abbc :aabc) (:abbc :abba) (:abba :aabc))) 
(fep-graph *pairs*)

||#


