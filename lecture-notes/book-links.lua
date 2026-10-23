-- Links from lecture notes into the book.
--
-- Write [link text](book:prp-rt-nogo) for any id the book defines ({#sec-x},
-- {#prp-x}, {#eq-x}, ...). The filter finds the book chapter that defines the id
-- and turns the target into
--   https://gustavdelius.github.io/affine_toda/<chapter path>.html#<id>
-- so the link survives a renumbering or a move of the chapter. Give the link a
-- descriptive text; never type a book number. An id that the book does not
-- define is left as it is, with a warning; check-refs.py, run on the notes
-- project, fails on it.
--
-- The book is the repository root: the nearest ancestor of the notes project
-- that holds code/. Book files are the .qmd files outside
-- lecture-notes/, briefs/, referee_scratch/ and build or tool directories, the
-- same rule as in check-refs.py.

local base = "https://gustavdelius.github.io/affine_toda/"
local skip = {
  ["lecture-notes"] = true, briefs = true, referee_scratch = true, _site = true,
  _book = true, [".quarto"] = true, [".git"] = true, [".claude"] = true,
  code = true, foundations_code = true, soliton_rmatrix_code = true,
  node_modules = true,
}

local function is_dir(path)
  local ok = pcall(pandoc.system.list_directory, path)
  return ok
end

local function project_dir()
  if quarto and quarto.project and quarto.project.directory then
    return quarto.project.directory
  end
  return pandoc.system.get_working_directory()
end

local function find_repo(dir)
  local d = dir
  for _ = 1, 20 do
    if is_dir(pandoc.path.join({ d, "code" })) then return d end
    local parent = pandoc.path.directory(d)
    if parent == d then break end
    d = parent
  end
  return nil
end

local repo = find_repo(project_dir())
local ids = nil

local function scan(dir, rel)
  for _, name in ipairs(pandoc.system.list_directory(dir)) do
    local path = pandoc.path.join({ dir, name })
    local relpath = rel == "" and name or (rel .. "/" .. name)
    if name:sub(-4) == ".qmd" then
      local f = io.open(path, "r")
      if f then
        local text = f:read("a")
        f:close()
        for id in text:gmatch("{#(%a+%-[%w_%-]+)") do
          if not ids[id] then ids[id] = relpath:sub(1, -5) .. ".html" end
        end
      end
    elseif not skip[name] and name:sub(1, 1) ~= "." and is_dir(path) then
      scan(path, relpath)
    end
  end
end

local function book_ids()
  if ids then return ids end
  ids = {}
  if repo then scan(repo, "") end
  return ids
end

local function Link(el)
  local id = el.target:match("^book:([%w_%-]+)$")
  if not id then return nil end
  local page = book_ids()[id]
  if not page then
    io.stderr:write("WARNING (book-links.lua): the book defines no id " .. id .. "\n")
    return nil
  end
  el.target = base .. page .. "#" .. id
  return el
end

return { { Link = Link } }
