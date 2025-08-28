function Qualite=psnr(ImRef,ImDis)
  Qualite=-10*log10(mean2((ImRef-ImDis).^2)+eps);
