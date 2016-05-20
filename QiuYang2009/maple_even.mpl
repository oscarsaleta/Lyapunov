fname:=cat("/home/osr/degree",taskId,"solved.txt");
res:=solve(taskArguments[1]=0,a);
print(res);
save res, fname;
