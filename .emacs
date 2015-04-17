;; .emacs
;; load viper
(setq viper-mode t)
(require 'viper)

(setq tex-dvi-view-command "xdvi")
(add-hook 'LaTeX-mode-hook (function turn-on-reftex))
(setq reftex-plug-into-AUCTeX t)

;; save place in file
(setq save-place-file "~/.emacs.d/saveplace") ;; keep my ~/ clean
(setq-default save-place t)                   ;; activate it for all buffers
(require 'saveplace)                          ;; get the package

;;;;;;;;;;;;;;;;;;;;;;;;; TABS
(add-to-list 'load-path "~/.emacs.d/elisp/")
(require 'smart-tabs-mode) 
(autoload 'smart-tabs-mode "smart-tabs-mode"
    "Intelligently indent with tabs, align with spaces!")
(autoload 'smart-tabs-mode-enable "smart-tabs-mode")
(autoload 'smart-tabs-advice "smart-tabs-mode")
(autoload 'smart-tabs-insinuate "smart-tabs-mode")
;;; this sets tab-width to 4
(setq-default indent-tabs-mode nil)
(setq-default tab-width 4)
(setq indent-line-function 'insert-tab)
;;; this insunates smart tabs -- must come after tabs-mode?
(smart-tabs-insinuate 'c 'c++ 'javascript )
;; highlighting (doesnt work?)
(setq-default highlight-tabs)
(setq-default highlight-trailing-whitespace)

;;;;change emacs indent
;;(setq c-default-style "cc-mode"
;;      c-basic-offset 4)
;;;;(setq split-height-threshold 0)
;;;;(setq split-width-threshold 0)
;;;;This will make spaces the indent character, and use 4 spaces per
;;;;indent level, for C, C++, and Objective C:
;;(setq c-mode-hook
;;    (function (lambda ()
;;                (setq c-indent-level 4))))
;;                (setq indent-tabs-mode nil)




;;; uncomment this line to disable loading of "default.el" at startup
;; (setq inhibit-default-init t)

;; turn on font-lock mode
(when (fboundp 'global-font-lock-mode)
  (global-font-lock-mode t))

;; enable visual feedback on selections
;(setq transient-mark-mode t)

;; default to better frame titles
(setq frame-title-format
      (concat  "%b - emacs@" (system-name)))

;; default to unified diffs
(setq diff-switches "-u")

;; always end a file with a newline
;(setq require-final-newline 'query)

(set-face-foreground 'font-lock-comment-face "Red" )
(set-variable font-lock-comment-face 'font-lock-comment-face) 


(custom-set-variables
  ;; custom-set-variables was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 '(inhibit-startup-screen t)
 '(vc-follow-symlinks t))
(custom-set-faces
  ;; custom-set-faces was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 '(sh-heredoc-face ((((class color) (background light)) (:foreground "green"))) t))
 ;;'(sh-heredoc-face ((((class color) (background light)) (:foreground "green")))))
(setq auto-mode-alist (cons '("\\.bib$" . text-mode) auto-mode-alist))


;;(add-hook 'c-mode-common-hook '(lambda ()
;;      (local-set-key (kbd "RET") 'newline-and-indent)))

(setq column-number-mode 1)

(setq sentence-end-double-space nil)

;;(hi-lock-exclude-modes )

		

(setq mouse-wheel-scroll-amount '(1 ((shift) . 1))) ;; one line at a time
;;(setq scroll-step 1)
(setq scroll-conservatively 10000)

;; These tell emacs to associate certain filename extensions with 
;; certain modes.  I use cc-mode.el (c++-mode) for C as well as C++
;; code.  It is fairly all-encompassing, also working with other C-like
;; languages, such as Objective C and Java. 
(setq auto-mode-alist (cons '("\\.m$" . octave-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.text$" . text-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.txt$" . text-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.doc$" . text-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.awk$" . awk-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.perl$" . perl-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.plx$" . perl-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.pl$" . perl-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.C$" . c-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.cc$" . c++-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.c$" . c-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.def$" . c++-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.h$" . c-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.cpp$" . c++-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.cxx$" . c++-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.tcl$" . tcl-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.sh$" . shell-script-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.cray$" . shell-script-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.zsh$" . shell-script-mode)   auto-mode-alist))
(setq auto-mode-alist (cons '("\\.pyx$" . python-mode)   auto-mode-alist))
(setq auto-mode-alist (cons '("\\.pxd$" . python-mode)   auto-mode-alist))

(add-hook 'sh-mode-hook
          (lambda ()
             (define-key sh-mode-map "\C-c\C-c" 'comment-region)))
			  
;;--------------------------------------------------------------------   
;; Lines enabling gnuplot-mode                                          

;; move the files gnuplot.el to someplace in your lisp load-path or    
;; use a line like                                                    
;;  (setq load-path (append (list "/path/to/gnuplot") load-path))    
(setq load-path (append (list "~/.emacs.d/gnuplot-mode") load-path))

;; these lines enable the use of gnuplot mode                       
  (autoload 'gnuplot-mode "gnuplot" "gnuplot major mode" t)
  (autoload 'gnuplot-make-buffer "gnuplot" "open a buffer in gnuplot mode" t)

;; this line automatically causes all files with the .gp extension to
;; be loaded into gnuplot mode                                      
  (setq auto-mode-alist (append '(("\\.gp$" . gnuplot-mode)) auto-mode-alist))
;;  (setq auto-mode-alist (append '(("\\.gnu$" . gnuplot-mode)) auto-mode-alist))

;; This line binds the function-9 key so that it opens a buffer into
;; gnuplot mode                                                    
  (global-set-key [(f9)] 'gnuplot-make-buffer)

;; end of line for gnuplot-mode              
;;--------------------------------------------------------------------


