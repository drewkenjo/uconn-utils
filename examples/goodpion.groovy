import org.jlab.io.hipo.HipoDataSource
import org.jlab.clas.physics.LorentzVector
import uconn.utils.pid.stefan.PionCandidate
import uconn.utils.pid.stefan.PionCandidate.Cut

def banknames = ['RUN::config', 'REC::Particle','REC::Traj']

args.each{fname->
  def reader = new HipoDataSource()
  reader.open(fname)

  while(reader.hasEvent()) {
    def event = reader.getNextEvent()
    if(banknames.every{event.hasBank(it)}) {
      def (runb,partb,trajb) = banknames.collect{event.getBank(it)}

      // all particles that pass pi+ cuts from PionCandidate class
      def pips0 = (0..<partb.rows()).collect{ipart->PionCandidate.getPionCandidate(ipart, partb, trajb)}.findAll{it.ispip()}
      pips0.each{candidate->
        println("good pi+ is found")
        def pip = candidate.getLorentzVector()
        println("pi+ momentum = "+pip.p())
      }

      //specify the list of pi+ cuts you want to apply
      def pips1 = (0..<partb.rows()).collect{ipart->PionCandidate.getPionCandidate(ipart, partb, trajb)}.findAll{it.ispip(Cut.PION_PID, Cut.DELTA_VZ)}
      pips1.each{candidate->
        println("good pi+ is found")
        def pip = candidate.getLorentzVector()
        println("pi+ momentum = "+pip.p())
      }

    }
  }

  reader.close()
}


