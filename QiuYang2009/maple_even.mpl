fname:=cat("degree",taskId,"solved.txt");
res:=solve(taskArguments[1]=0,a);
save res, fname;
