-- Put each part's introduction on the part's title page in the PDF.
--
-- Quarto turns a part's title into a paragraph "\part{...}" and leaves the
-- text of part.qmd after it, so KOMA-Script sets that text on the page after
-- the part page. Moving the text into \setpartpreamble, placed before \part,
-- sets it under the part title instead. Runs at post-render, after Quarto has
-- created the \part paragraphs.

local function starts_with_raw_latex(block)
  return block.t == "Para" and #block.content > 0
    and block.content[1].t == "RawInline"
end

local function is_part(block)
  return starts_with_raw_latex(block)
    and block.content[1].text:match("^\\part{") ~= nil
end

local intro_types = { Para = true, Plain = true, BulletList = true,
                      OrderedList = true, BlockQuote = true }

-- Applied to every list of blocks, since Quarto wraps each part in a Div.
local function move_intros(blocks)
  local out, i = pandoc.Blocks({}), 1
  while i <= #blocks do
    local block = blocks[i]
    out:insert(block)
    i = i + 1
    if is_part(block) then
      local intro = pandoc.Blocks({})
      while i <= #blocks and intro_types[blocks[i].t]
          and not starts_with_raw_latex(blocks[i]) do
        intro:insert(blocks[i])
        i = i + 1
      end
      if #intro > 0 then
        local tex = pandoc.write(pandoc.Pandoc(intro), "latex")
        out:insert(#out, pandoc.RawBlock("latex",
          "\\setpartpreamble[u][0.8\\textwidth]{\\vspace{2em}\n" .. tex .. "\n}"))
      end
    end
  end
  return out
end

if FORMAT:match("latex") then
  return { { Blocks = move_intros } }
end
