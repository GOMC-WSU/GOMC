  //std::cout << newMolPos << std::endl;
  //std::cout << newCOM << std::endl;
    std::cout << "B4" << std::endl;
  XYZArray testNMP = newMolPos;
  XYZ testCOM = newCOM;
    std::cout << "B4 POS" << std::endl;

  for (int i = 0; i < testNMP.Count(); ++i){
    std::cout << newMolPos.x[i] << ", " << newMolPos.y[i]
     << ", " << newMolPos.z[i] << std::endl;
  }

  for (int i = 0; i < testNMP.Count(); ++i){
    std::cout << coordCurrRef.x[pStart+i] << ", " << coordCurrRef.y[pStart+i]
     << ", " << coordCurrRef.z[pStart+i] << std::endl;
  }
  std::cout << "B4 COM" << std::endl;

    std::cout << newCOM.x << ", " << newCOM.y
     << ", " << newCOM.z << std::endl;
  
  coordCurrRef.TranslateRand(testNMP, testCOM, pStart, pLen,
                             m, b, moveSetRef.Scale(b, mv::DISPLACE, mk));
    std::cout << "TEST POS" << std::endl;

  for (int i = 0; i < testNMP.Count(); ++i){
    std::cout << testNMP.x[i] << ", " << testNMP.y[i]
     << ", " << testNMP.z[i] << std::endl;
  }
  std::cout << "TEST COM" << std::endl;

    std::cout << testCOM.x << ", " << testCOM.y
     << ", " << testCOM.z << std::endl;

  std::cout << "AFTER POS" << std::endl;

  for (int i = 0; i < testNMP.Count(); ++i){
    std::cout << newMolPos.x[i] << ", " << newMolPos.y[i]
     << ", " << newMolPos.z[i] << std::endl;
  }
  std::cout << "AFTER COM" << std::endl;

    std::cout << newCOM.x << ", " << newCOM.y
     << ", " << newCOM.z << std::endl;
  

  exit(1);