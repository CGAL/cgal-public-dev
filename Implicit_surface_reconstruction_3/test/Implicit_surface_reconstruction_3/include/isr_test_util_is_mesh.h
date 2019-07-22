bool is_mesh(std::string input_filename)
{
  std::string extension = input_filename.substr(input_filename.find_last_of('.'));
  if (extension == ".off" || extension == ".OFF" || extension == ".stl" || extension == ".STL" || extension == ".obj" || extension == ".OBJ") 
  {
    return true;
  }
  return false;
}