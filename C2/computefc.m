function fc = comptuefc(concrete, ec, alternative)
  % This function computes the concrete's stress in function of the concret's deformation using
  % concrete: struct with the following concrete's properties:
  %    E: Young Modulus
  %    eo: 
  % ec: deformation of the fiber (scalar or vector)
  % alternative: char with the name of the method to estimate the stress ('mander', '...')

  % Choose alternative
  if isequal(lower(alternative), 'mander')
    
  elseif isequal(lower(alternative), '')
    
  else

  end
end

