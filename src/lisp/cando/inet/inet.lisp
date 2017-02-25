(in-package :inet)


(defun download-pdb (accession-code pathname)
  "Download the PDB file from RCSB and put it in the (pathname)"
  (let ((url (concatenate 'string "http://files.rcsb.org/download/" accession-code ".pdb1")))
    (trivial-http:http-download url pathname)))
