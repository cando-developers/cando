(in-package :inet)


(defun download-pdb (accession-code &optional pathname)
  "Download the PDB file from RCSB and put it in the (pathname)"
  (let ((accession-name (string accession-code)))
    (unless pathname (setf pathname (merge-pathnames (pathname (concatenate 'string accession-name ".pdb")))))
    (let ((url (concatenate 'string "http://files.rcsb.org/download/" accession-name ".pdb1")))
      (trivial-http:http-download url pathname))))
