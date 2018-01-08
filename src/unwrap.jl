"""

unwrapped = unwrap(phase,inplace=false)

Given an input phase `phase` (an array), returns an unwrapped array.
"""

function unwrap(phase, inplace=false)

  unwrapped = inplace ? phase : copy(phase)
  for i in 2:length(phase)
    while unwrapped[i] - unwrapped[i-1] >= π
      unwrapped[i] -= 2π
    end
    while unwrapped[i] - unwrapped[i-1] <= -π
      unwrapped[i] += 2π
    end
  end
  return unwrapped
end

unwrap!(phase) = unwrap(phase, true)
