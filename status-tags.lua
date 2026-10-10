-- Website only: turn the bold status labels, such as "**[Theorem]**" or
-- "**[Proposition, sketch]**", into colour-coded tags.
--
-- A label is a Strong in brackets that opens a paragraph or list item and
-- starts with one of the status words below; it appears in theorem-like
-- blocks, in lists of results, and in the legend in the introduction. The
-- brackets are dropped and the label is wrapped in a Span with classes
-- "status-tag status-<kind>", which theme.scss styles. Runs at pre-ast, before
-- Quarto restructures theorem blocks. The PDF keeps the bracketed labels.

local status_words = { Theorem = true, Established = true, Proposition = true,
                       Numerical = true, Conjecture = true, Part = true }

local function kind_of(label)
  if label:match("^Part") then
    return "mixed"
  elseif label:find("sketch") or label:match("^Numerical") then
    return "partial"
  elseif label:match("^Theorem") or label:match("^Established") then
    return "proved"
  elseif label:match("^Conjecture") then
    return "conjecture"
  end
  return "mixed"
end

local function tag_label(block)
  local strong = block.content[1]
  if not strong or strong.t ~= "Strong" then
    return nil
  end
  local inlines = pandoc.Inlines(strong.content)
  local n = #inlines
  local first, last = inlines[1], inlines[n]
  if not (first and first.t == "Str" and first.text:sub(1, 1) == "["
          and last.t == "Str" and last.text:sub(-1) == "]") then
    return nil
  end
  if n == 1 then
    inlines[1] = pandoc.Str(first.text:sub(2, -2))
  else
    inlines[1] = pandoc.Str(first.text:sub(2))
    inlines[n] = pandoc.Str(last.text:sub(1, -2))
  end
  local label = pandoc.utils.stringify(inlines)
  if not status_words[label:match("^%a+") or ""] then
    return nil
  end
  block.content[1] = pandoc.Span(inlines,
    { class = "status-tag status-" .. kind_of(label) })
  return block
end

if not quarto.doc.is_format("html") then
  return {}
end

return { { Para = tag_label, Plain = tag_label } }
