-- Print cross-references to the hypothesis floats (ids hyp-D, hyp-X1, ...)
-- by their letter, as "Hypothesis D", instead of Quarto's per-chapter number.
-- The link keeps the anchor of the float, so the reference still resolves.
-- check-refs.py still checks that every @hyp- id is defined.
function Cite(el)
  if #el.citations == 1 and el.citations[1].id:match('^hyp%-') then
    local name = el.citations[1].id:match('^hyp%-(.+)$')
    return pandoc.Link({pandoc.Str('Hypothesis ' .. name)}, '#' .. el.citations[1].id)
  end
end
