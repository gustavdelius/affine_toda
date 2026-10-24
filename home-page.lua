-- Keep the home page of the site out of the PDF.
--
-- index.qmd is the home page: it says what the site holds and compares the
-- lecture courses, inside a div shown in HTML only. Quarto still takes the
-- page's first heading as the title of a book chapter, so without this filter
-- the PDF would open with an empty chapter "About this site" before the
-- preface. The filter removes that heading from the LaTeX output.

function Header(el)
  if el.identifier == "sec-about-site" and quarto.doc.is_format("latex") then
    return {}
  end
end
