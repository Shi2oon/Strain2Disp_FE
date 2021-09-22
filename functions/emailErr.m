function emailErr(err,Dirunique,email_add)
% write the error to string for email % first line: message
message{1}=sprintf('for data %s %s\n',Dirunique,err.message);

% following lines: stack
  for e=1:length(err.stack)
     message{e+1}=sprintf('%s in %s at %i\n',err.stack(e).name,...
							err.stack(e).file,err.stack(e).line);
  end
  
gmail(email_add,['Error!: xEBSD2ABAQUS on ' getenv('COMPUTERNAME') ' for ' ...
		getenv('USERNAME')],message)
rethrow(err);
end