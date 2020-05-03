function exportAsPdf(handle, name_figure)

set(handle,'Units','inches');
screenposition = get(handle,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',[screenposition(3:4)]);
print (name_figure, '-dpdf', '-painters')
set(gcf,'renderer','Painters')
system(['pdfcrop ',name_figure,'.pdf ',name_figure,'.pdf']);

end