function fingerTappingTest1()
  clear all
  clc
  [file,dir] = uigetfile('.txt','Importar o arquivo',"multiselect","on");
  [m,n]=size(file);

  for z = 1:n
    z

    arquivo = [dir file{z}];
    info = importdata(arquivo);
    resolution1 = info.textdata{7,1};
    resolution2 = info.textdata{8,1};
    width = str2num(resolution1(1,8:end))*0.026458;
    height = str2num(resolution2(1,9:end))*0.026458;
    c = csvread(arquivo);
    t = c(10:end,1)/1000;
    x = c(10:end,2)*0.026458;
    y = c(10:end,3)*0.026458;
    f = c(10:end,4);
    v = sum(sqrt(x.^2+y.^2));

    fim = ceil(t(end,1));
    tempo = 1:1:fim;
    c = 0;
    for i = 1:length(tempo)
      for j = 1:length(t)
        if t(j,1) < tempo(1,i)
          density(1,i) = c;
        endif
      endfor
      c = 0;
    endfor

    for i = 1:length(t)-1
      freq(i,1) = 1/(t(i+1,1)-t(i,1));
    endfor
    n_errors = 0;

    g1 = 1;
    g2 = 1;

    for i = 1:length(t)
      if f(i,1) == 1
        group_1(g1,1:2) = [x(i,1) y(i,1)];
        g1 = g1 + 1;
      elseif f(i,1) == 2
        group_2(g2,1:2) = [x(i,1) y(i,1)];
        g2 = g2 + 1;
      endif
    endfor

    if g1>1 && g2>1
      [a1,b1,D1,d1] = ellipseModel(group_1(:,1),group_1(:,2));
      [a2,b2,D2,d2] = ellipseModel(group_2(:,1),group_2(:,2));
    else
      [a1,b1,D1,d1] = ellipseModel(group_1(:,1),group_1(:,2));
    endif

    if f(1,1) == 0
      n_errors = n_errors + 1;
    else
    endif

    for i = 2:length(x)

      if g1>1 && g2==1
        if f(i,1) == 1
        else
          n_errors = n_errors + 1;
        endif
      else
        if f(i,1) == 1 && f(i-1,1) ~= 1
        elseif f(i,1) == 2 && f(i-1,1) ~= 2
        else
          n_errors = n_errors + 1;
        endif
      endif

    endfor

    if g1>1 && g2>1
      [vec,val]=eig(cov(group_1(:,1),group_1(:,2)));
      Area_1 = pi*prod(2.4478*sqrt(svd(val)));
      [vec,val]=eig(cov(group_2(:,1),group_2(:,2)));
      Area_2 = pi*prod(2.4478*sqrt(svd(val)));

      resultado = [n_errors;max(freq);min(freq);median(freq);sum(v);Area_1;Area_2;abs((Area_2)-(Area_1));length(x);length(x)/max(t/1000);D1;d1;pi*(D1)*(d1);D2;d2;pi*(D2)*(d2)];
    else

      [vec,val]=eig(cov(group_1(:,1),group_1(:,2)));
      Area_1 = pi*prod(2.4478*sqrt(svd(val)));

      resultado = [n_errors;max(freq);min(freq);median(freq);sum(v);Area_1;length(x);length(x)/max(t/1000);D1;d1;pi*(D1)*(d1)];
    endif

    dir2 = dir;
    file{z} = strrep(file{z},".txt"," resultados.txt");
    dir2 = strrep(dir2,"testes","resultados");
    arquivo = [dir2 file{z}]
    dlmwrite(arquivo,resultado,"\t");

    keepvars = {'file', 'dir', 'n'};
    clearvars('-except', keepvars{:});

  endfor
  beep();
