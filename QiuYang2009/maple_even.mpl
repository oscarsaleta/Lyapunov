fname:=cat("degree",taskId,"solved.txt");
res:=solve(taskArgs[1]=0,a);
print(res);
save res, fname;
