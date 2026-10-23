-- Link the scripts that the text cites to their source on GitHub.
--
-- An inline code span such as `code/boundary_coideal.py`,
-- `code/hirota/fold.py --grid 768` or `code/hirota/` becomes
-- a link when its first whitespace-delimited token names an existing file (or,
-- for a token ending in "/", an existing directory) in the repository. The token
-- is tried relative to the repository root (the project root, or for a nested
-- project such as lecture-notes/defects/ the enclosing repository), then
-- relative to code/, so
-- that a table may use bare names such as `defect_axioms.py`. Any other code span is left alone, as is code that
-- already sits inside a link. Works for HTML and PDF.

local base = "https://github.com/gustavdelius/affine_toda/"

local function find_root()
  if quarto and quarto.project and quarto.project.directory then
    return quarto.project.directory
  end
  local input = PANDOC_STATE and PANDOC_STATE.input_files and PANDOC_STATE.input_files[1]
  if not input or input == "" then return nil end
  if input:sub(1, 1) ~= "/" then
    input = pandoc.path.join({ pandoc.system.get_working_directory(), input })
  end
  local dir = pandoc.path.directory(input)
  for _ = 1, 20 do
    local f = io.open(pandoc.path.join({ dir, "_quarto.yml" }), "r")
    if f then f:close(); return dir end
    local parent = pandoc.path.directory(dir)
    if parent == dir then break end
    dir = parent
  end
  return nil
end

local cache = {}

-- "file", "dir" or false
local function kind(path)
  local hit = cache[path]
  if hit ~= nil then return hit end
  local result = false
  local f = io.open(path, "r")
  if f then
    local _, err, code = f:read(0)
    if err and code == 21 then
      result = "dir"            -- EISDIR on Linux and macOS
    else
      result = "file"
    end
    f:close()
  end
  cache[path] = result
  return result
end

-- The repository root: the project directory, or, for a project nested in the
-- repository such as lecture-notes/defects/, the nearest ancestor that holds
-- code/. Paths are resolved and linked relative to it.
local function find_repo(dir)
  if not dir then return nil end
  local d = dir
  for _ = 1, 20 do
    if kind(pandoc.path.join({ d, "code" })) == "dir" then return d end
    local parent = pandoc.path.directory(d)
    if parent == d then break end
    d = parent
  end
  return dir
end

local root = find_repo(find_root())

local function resolve(token)
  if not root or token == "" or token:sub(1, 1) == "/" or token:find("..", 1, true) then
    return nil
  end
  local want_dir = token:sub(-1) == "/"
  local rel = want_dir and token:sub(1, -2) or token
  if rel == "" then return nil end
  for _, prefix in ipairs({ "", "code/" }) do
    local k = kind(root .. "/" .. prefix .. rel)
    if k == (want_dir and "dir" or "file") then
      return prefix .. rel, k
    end
  end
  return nil
end

local function Code(el)
  local token = el.text:match("^%s*(%S+)")
  if not token then return nil end
  local rel, k = resolve(token)
  if not rel then return nil end
  local url = base .. (k == "dir" and "tree/main/" or "blob/main/") .. rel
  return pandoc.Link({ el }, url), false
end

-- Do not look inside existing links.
local function Link(el)
  return el, false
end

return { { traverse = "topdown", Link = Link, Code = Code } }
