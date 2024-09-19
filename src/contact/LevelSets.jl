function plane(x, y)
  return y - x[2, :]
end

function corner(x, x_loc, y_loc)
  return min(x[1, :] - x_loc, x[2, :] - y_loc)
end

function sphere(x, x_loc, y_loc, R)
  r = sqrt((x[1, :] - x_loc)^2 + (x[2, :] - y_loc)^2) - R
  return r
end

function combined(x, ls1, ls2)
  return min(ls1(x), ls2(x))
end
