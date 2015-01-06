#! define assertions that can be disabled if performance becomes an issue...
macro ndebug()
  return false;
end

macro dassert(ex, msg)
  return quote
    if !(@ndebug)
      @assert(ex, msg);
    end
  end
end