import org.jlab.io.hipo.HipoDataSource
import org.jlab.clas.physics.LorentzVector
import uconn.utils.pid.stefan.ElectronCandidate
import uconn.utils.pid.stefan.ElectronCandidate.Cut

def banknames = ['RUN::config', 'REC::Particle','REC::Calorimeter','REC::Cherenkov','REC::Traj']

args.each{fname->
  def reader = new HipoDataSource()
  reader.open(fname)

  while(reader.hasEvent()) {
    def event = reader.getNextEvent()
    if(banknames.every{event.hasBank(it)}) {
      def (runb,partb,calb,ccb,trajb) = banknames.collect{event.getBank(it)}

      def ipart = 0
      def candidate = ElectronCandidate.getElectronCandidate(ipart, partb, calb, ccb, trajb)

      //test for all cuts specified in ElectronCandidate class
      if(candidate.iselectron()) {
        println("good electron is found")
        def ele = candidate.getLorentzVector()
        println("electron momentum = "+ele.p())
      }

      //specify the list of cuts you want to apply
      if(candidate.iselectron(Cut.ELE_PID, Cut.CC_NPHE, Cut.DC_VERTEX)) {
        println("good electron is found")
        def ele = candidate.getLorentzVector()
        println("electron momentum = "+ele.p())
      }

    }
  }

  reader.close()
}


