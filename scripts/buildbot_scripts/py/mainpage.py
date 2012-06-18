#================================== NEW ================================

from buildbot.status.web.base import HtmlResource
import mainwork
import os

class mainPage(HtmlResource):
    title = "Define Experiment"
    test_Info=""
    change_Info=""
    exp_Info=""

    def __init__(self, categories=None):
       HtmlResource.__init__(self)
#       self.categories = categories
       self.putChild("work", mainwork.MainWork())

    def main_Change(self, request):
      global change_Info
      data =  '<meta HTTP-EQUIV="REFRESH" content="0; url=../">'
      
      if change_Info == "new":

# Create a emty Experiment List type

        ExpList = []
        ExpListInfo = []

# Read all files in the "icon_run" directory and write the filename which starts with "exp.test_" into the
# Experiment List type

        for filename in os.listdir("icon_run"):
	  if filename.find("exp.test_") == 0:
            ExpDict = {"Description" : " ","Model" : " ","Grid" : " "}
            ExpDict["ExpName"] = filename.lstrip("exp.")
            
            file = open("icon_run/" + filename)
            while 1:
              line = file.readline()
              if not line:
                break
#              print "INFO: " + line
              if line.find("_Description_") >= 0:
		tmp,Info = line.split("_Description_",1)
                Info.replace("\n","")
                ExpDict["Description"] = Info.strip(" ")
                
              if line.find("_Model_") >= 0:
		tmp,Info = line.split("_Model_",1)
                Info.replace("\n","")
                ExpDict["Model"] = Info.strip(" ")
                
              if line.find("_Grid_") >= 0:
		tmp,Info = line.split("_Grid_",1)
                Info.replace("\n","")
                ExpDict["Grid"] = Info.strip(" ")
              
            file.close
            ExpList.append(filename.lstrip("exp."))
            ExpListInfo.append(ExpDict)
            
            
# Sort the Experiment List type

        ExpList.sort()
	
        data  = '''
<h1>ICON Buildbot</h1>
Buildbot starts every night a test suite on a set of builders differing in machines, compilers 
and parallelization setup. A standard Buildbot test unpacks the model code and scripts from the  
repositors, configures and compiles the codes for the selected builder, and runs a series of 
standard test experiments and related post-processings. Figures resulting from the 
post-processings are archived and are uploaded to web pages for comparison with reference 
figures.<p>
	'''
	
	data += '''
<h2>Status results</h2>
<table border="1" cellspacing="0" cellpadding="2">
   <thead>
     <tr valign="baseline">
       <td><i>Display</i></td>
       <td><i>What it shows</i></td>
     </tr>
   </thead>
   <tfoot></tfoot>
   <tbody>
     <tr valign="baseline">
       <td><a href="waterfall?show_events=false&show_time=604800">Waterfall</a></td>
       <td>Status of each step of the Buildbot tests for the past 7 days</td>
     </tr>

     <tr valign="baseline">
       <td><a href="grid?width=10">Grid</a></td>
       <td>Status of the complete Builbot test for the 10 latest tested revisions</td>
     </tr>

     <tr valign="baseline">
       <td><a href="one_box_per_builder">Latest build</a></td>
       <td>Status of the complete Builbot test for the latest test on each builder</td>
     </tr>

     <tr valign="baseline">
       <td><a href="one_line_per_build">List</a></td>
       <td>Status of the complete Builbot test for the latest 15 tests</td>
     </tr>

     <tr valign="baseline">
       <td><a href="buildslaves">Build slaves</a></td>
       <td>Information on machines ("slaves"), builders, etc.</td>
     </tr>
   </tbody>
   </table>
        '''

	data += '''
<h2>Eperiment results</h2>
<table border="1" cellspacing="0" cellpadding="2">
  <thead>
    <tr valign="baseline">
      <td width="300"><i>Experiment (exp)</i></td>
      <td width="380"><i>Description</i></td>
      <td width="200"><i>Model</i></td>
      <td width="120"><i>Grid</i></td>
    </tr>
  </thead>
  <tfoot>
  </tfoot>
  <tbody>
        '''
        
# Build HTML Table info

        for exp in ExpList:
          for e in ExpListInfo:
            if exp == e.get('ExpName'):
	      break
	      
 	  data += '''
<tr valign="baseline">
          '''
          data += "<td><a href=\"plot?exp=" + e.get('ExpName') + "&modus=nightly\">" + e.get('ExpName') + "</a></td>"
          data += "<td>" + e.get('Description') + "</td>"
          data += "<td>" + e.get('Model') + "</td>"
          data += "<td>" + e.get('Grid') + "</td>"
          data += "</tr>"
            

	data += '''
 	  </tbody>
	  </table>
	'''

        data += '''
          <h3>About Buildbot</h3>
	  <ul>
	    <li><a href="about">Buildbot version used here</a></li>
	    <li>Buildbot home page: <a href="http://trac.buildbot.net" target="_blank">http://trac.buildbot.net</a></li>
	  </ul>
	'''

      return data
  
    def plot_exp_info(self, request):
      global test_Info
      global change_Info
      global exp_Info
      data = 'Page exp test'+"<br>"
      data += "test_Info= "+test_Info+"<br>"
      data += "change_Info= "+change_Info+"<br>"
      data += "exp_Info= "+exp_Info+"<br>"
      return data
      
    def plot_exp(self, request):
      global test_Info
      global change_Info
      global exp_Info
      data = 'Page exp'+"<br>"
      data += "test_Info= "+test_Info+"<br>"
      data += "change_Info= "+change_Info+"<br>"
      data += "exp_Info= "+exp_Info+"<br>"
      return data
    
    def get_info(self, request):
        global test_Info
        global change_Info
        global exp_Info
	
        test_Info="NotSet"
        change_Info="NotSet"
        exp_Info="NotSet"
        
	if "status" in request.args:
          try:
            change_Info = request.args["status"][0]
          except ValueError:
            pass
	
        if "exp" in request.args:
          try:
            exp_Info = request.args["exp"][0]
          except ValueError:
            pass
	
        if "test" in request.args:
          try:
            test_Info = request.args["test"][0]
          except ValueError:
            pass
	  	
        return None
    
    def head(self, request):
        return ""
    
    #def head(self, request):
        #head = '<meta HTTP-EQUIV="REFRESH" content="0; url=../">'
        #return head
    
    def body(self, request):
        global test_Info
        global change_Info
        global exp_Info
	self.get_info(request)

	if exp_Info != "NotSet":
	  if test_Info != "NotSet":
            data = self.plot_exp_test(request)
	  else:
	    data = self.plot_exp(request) 
	else:
	  data = self.main_Change(request)	
	  
        return data
#================================== NEW ================================
