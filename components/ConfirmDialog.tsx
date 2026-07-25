import React, { useEffect, useRef } from 'react';

/**
 * ConfirmDialog — themed confirmation built on the native `<dialog>` element.
 *
 * `showModal()` gives us focus trapping, inertness of the page behind, Esc to
 * dismiss, and the top-layer stacking that would otherwise need a portal and a
 * hand-rolled focus trap. We take all of that and restyle only the paint.
 *
 * `window.confirm()` was the obvious alternative and is the wrong one: it
 * renders OS chrome in the middle of a near-black app, which is the same
 * foreign-control problem the themed Select and scrollbars exist to fix.
 */
interface ConfirmDialogProps {
  open: boolean;
  title: string;
  body: string;
  confirmLabel: string;
  cancelLabel?: string;
  /** Destructive actions get the red treatment on the primary button. */
  destructive?: boolean;
  onConfirm: () => void;
  onCancel: () => void;
}

const ConfirmDialog: React.FC<ConfirmDialogProps> = ({
  open, title, body, confirmLabel, cancelLabel = 'Cancel',
  destructive = false, onConfirm, onCancel,
}) => {
  const ref = useRef<HTMLDialogElement>(null);

  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    if (open && !el.open) el.showModal();
    else if (!open && el.open) el.close();
  }, [open]);

  // Esc and backdrop dismissal both surface as `cancel`/`close`; route them
  // through onCancel so the caller's pending state is always cleared.
  useEffect(() => {
    const el = ref.current;
    if (!el) return;
    const onCancelEvent = (e: Event) => { e.preventDefault(); onCancel(); };
    el.addEventListener('cancel', onCancelEvent);
    return () => { el.removeEventListener('cancel', onCancelEvent); };
  }, [onCancel]);

  return (
    <dialog
      ref={ref}
      aria-labelledby="confirm-title"
      className="rg-dialog max-w-sm border border-edge bg-surface p-0 text-ink shadow-2xl"
      onClick={(e) => { if (e.target === ref.current) onCancel(); }}
    >
      <div className="p-4">
        <h2 id="confirm-title" className="text-sm font-semibold text-ink-strong">{title}</h2>
        <p className="mt-2 text-xs leading-relaxed text-ink-muted">{body}</p>
        <div className="mt-4 flex justify-end gap-2">
          <button
            onClick={onCancel}
            className="border border-edge bg-surface-raised px-3 py-1.5 text-xs text-ink-soft transition-colors hover:bg-surface-hover hover:text-ink-strong focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-brand/70"
          >
            {cancelLabel}
          </button>
          <button
            onClick={onConfirm}
            autoFocus
            className={`border px-3 py-1.5 text-xs font-medium text-ink-strong transition-colors focus-visible:outline-none focus-visible:ring-2 ${
              destructive
                ? 'border-danger bg-red-600 hover:bg-danger focus-visible:ring-danger-soft/70'
                : 'border-brand bg-brand-strong hover:bg-brand focus-visible:ring-brand-soft/70'
            }`}
          >
            {confirmLabel}
          </button>
        </div>
      </div>
    </dialog>
  );
};

export default ConfirmDialog;
