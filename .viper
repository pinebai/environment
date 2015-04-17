(setq viper-inhibit-startup-message 't)
(setq viper-expert-level '3)
(defadvice viper-maybe-checkout (around viper-svn-checkin-fix activate)
  "Advise viper-maybe-checkout to ignore svn files."
  (let ((file (expand-file-name (buffer-file-name buf))))
    (when (and (featurep 'vc-hooks)
	       (not (memq (vc-backend file) '(nil SVN))))
      ad-do-it)))
